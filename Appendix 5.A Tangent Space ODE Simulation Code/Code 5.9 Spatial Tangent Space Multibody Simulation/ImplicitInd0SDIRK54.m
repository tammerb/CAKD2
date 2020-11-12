function [v,vd,vdd,Lam,Jiter,R1n,J,JCond,h,nch]=...
    ImplicitInd0SDIRK54(n,tn,Vv,Vvd,Vvdd,LLam,Uu,q0,U,V,B,par,...
    h,hmax,nch,npar,SJDT,SMDT,STSDAT)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

%Runge-Kutta Coefficients
a=1/4;
aa=0.0625;
a1=[0.25,0,0,0,0;0.5,0.25,0,0,0;0.34,-0.04,0.25,0,0;0.272794117647059,...
    -0.0503676470588235,0.0275735294117647,0.25,0;1.04166666666667,...
    -1.02083333333333,7.8125,-7.08333333333333,0.25];
A1=a1*a1;
b=[1.041666666666667,-1.020833333333333,7.812500000000000,...
    -7.083333333333333,0.25];
B1=b*a1;
c=[1/4;3/4;11/20;1/2;1];
d1=[59/48,-17/96,225/32,-85/12,0];
D1=d1*a1;

Err=2;      %%Criteria for accepting time step
while Err>1

% Integration Jacobian Evaluation
tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
Lamnm=LLam(:,n-1);
u=Uu(:,n-1);

[Rvdd,Rvd,Rv,Rhs,qnm,Dnm]=JacobInd0(tnm,vnm,vdnm,vddnm,u,par,q0,V,U,B,...
    SJDT,SMDT,STSDAT);
Phiqnm=PhiqEval(tnm,qnm,SJDT,par);
J=[Rvdd+h*a*Rvd+(h^2)*aa*Rv,Phiqnm'];
JCond=cond(J);

%Third derivative calculation, with no explicit time dependence
M=MEval(qnm,SMDT,par);
Phiqnm=PhiqEval(tnm,qnm,SJDT,par);
EE=[M*Dnm,Phiqnm'];
xx=EE\Rhs;
Pvdd=[eye(nv),zeros(nv,nc)];
PLam=[zeros(nc,nv),eye(nc)];
vdddnm=Pvdd*xx;
Lamdnm=PLam*xx;

%Solution Estimates
vdd=vddnm+h*vdddnm; 
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
Lam=Lamnm+h*Lamdnm;

% Start Integration Step

% Solve Discretized Equations

% Stage 1
i=1;        %Set iteration counter
err=intol+1;

k1=vddnm;
z1dk=vdnm;
z1k=vnm+h*c(1)*vdnm;
while err>intol    
% Residual Calculation
v=z1k+(h^2)*aa*k1;
vd=z1dk+h*a*k1;

R1=ResidInd0(tn,v,vd,k1,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R11Norm=norm(R1);
end

% Quasi Newton Correction
x=-J\R1;
k1=k1+Pvdd*x;
Lam=Lam+PLam*x;
err=norm(R1);
i=i+1;

end
j1=i-1;

% Stage 2
i=1;        %Set solution iteration counter
err=intol+1;
k2=k1;
z2dk=vdnm+h*a1(2,1)*k1;
z2k=vnm+h*c(2)*vdnm+(h^2)*A1(2,1)*k1;
while err>intol   
% Residual Calculation
v=z2k+(h^2)*aa*k2;
vd=z2dk+h*a*k2;

R2=ResidInd0(tn,v,vd,k2,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R21Norm=norm(R2);
end

% Quasi Newton Correction
x=-J\R2;
k2=k2+Pvdd*x;
Lam=Lam+PLam*x;
err=norm(R2);
i=i+1;

end
j2=i-1;

% Stage 3
i=1;        %Set solution iteration counter
err=intol+1;
k3=k2;
z3dk=vdnm+h*(a1(3,1)*k1+a1(3,2)*k2);
z3k=vnm+h*c(3)*vdnm+(h^2)*(A1(3,1)*k1+A1(3,2)*k2);
while err>intol    
% Residual Calculation
v=z3k+(h^2)*aa*k3;
vd=z3dk+h*a*k3;

R3=ResidInd0(tn,v,vd,k3,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R31Norm=norm(R2);
end

% Quasi Newton Correction
x=-J\R3;
k3=k3+Pvdd*x;
Lam=Lam+PLam*x;
err=norm(R3);
i=i+1;

end
j3=i-1;

% Stage 4
i=1;        %Set solution iteration counter
err=intol+1;
k4=k3;
z4dk=vdnm+h*(a1(4,1)*k1+a1(4,2)*k2+a1(4,3)*k3);
z4k=vnm+h*c(4)*vdnm+(h^2)*(A1(4,1)*k1+A1(4,2)*k2+A1(4,3)*k3);
while err>intol    
% Residual Calculation
v=z4k+(h^2)*aa*k4;
vd=z4dk+h*a*k4;

R4=ResidInd0(tn,v,vd,k4,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R41Norm=norm(R2);
end

% Quasi Newton Correction
x=-J\R4;
k4=k4+Pvdd*x;
Lam=Lam+PLam*x;
err=norm(R4);
i=i+1;
end
j4=i-1;

% Stage 5
i=1;        %Set solution iteration counter
err=intol+1;
k5=k4;
z5dk=vdnm+h*(a1(5,1)*k1+a1(5,2)*k2+a1(5,3)*k3+a1(5,4)*k4);
z5k=vnm+h*c(5)*vdnm+...
    (h^2)*(A1(5,1)*k1+A1(5,2)*k2+A1(5,3)*k3+A1(5,4)*k4);
while err>intol    
% Residual Calculation
v=z5k+(h^2)*aa*k5;
vd=z5dk+h*a*k5;

R5=ResidInd0(tn,v,vd,k5,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R51Norm=norm(R2);
end

% Quasi Newton Correction
x=-J\R5;
k5=k5+Pvdd*x;
Lam=Lam+PLam*x;
err=norm(R5);
i=i+1;
end
j5=i-1;

%Calculate Solution

R1=[R11Norm;R21Norm;R31Norm;R41Norm;R51Norm];
R1n=norm(R1,inf);
jiter=[j1;j2;j3;j4;j5];
Jiter=norm(jiter,inf);

k=[k1,k2,k3,k4,k5];

%Evaluate v, vd, vdd, and Lam
v=vnm+h*vdnm+(h^2)*k*B1';
vd=vdnm+h*k*b';

[vdd,Lam,ECond]=Ind0ODEfunct(tn,v,vd,SMDT,STSDAT,SJDT,...
    u,q0,V,U,B,par);

if hvar==1  %variable step

%Evaluate imbedded solution for error estimate, vnhat and vdnhat
vnhat=vnm+h*vdnm+(h^2)*k*D1';
vdnhat=vdnm+h*k*d1';

%Evaluate  error with sciv=Atol*(1+abs(v(i))), scivd=Atol*(1+abs(vd(i)))...
%and p=5
vdiff=v-vnhat;
vddiff=vd-vdnhat;
Er=0;
i=1;
while i<=ngc-nc
sciv=Atol*(1+abs(v(i)));
scivd=Atol*(1+abs(vd(i)));
Er=Er+(vdiff(i)/sciv)^2+(vddiff(i)/scivd)^2;
i=i+1;
end
Err=sqrt(Er/(2*(ngc-nc)));

%Change step size

hopt=h*(1/Err)^(1/6);

if hopt<h
h=0.9*hopt;
nch=n;
if h<10^-5
h=10^-5;
Err=0.5;    
end
end

if hopt>h
    if n>nch+5
h=min([2*h;0.9*hopt]);
nch=n;
    end
end

if h>hmax
    h=hmax; 
end

end

if hvar==2
Err=0.5;
end

end

end



