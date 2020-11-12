function [v,vd,vdd,Lam,R1n,Jiter,JCond,h,nch]=...
    ImplicitInd0Trap(n,tn,npar,Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,h,hmax,nch,...
    SMDT,STSDAT,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

tnm=tn-h;
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
Lamnm=LLam(:,n-1);
unm=Uu(:,n-1);

Err=2;      %%Criteria for accepting time step
while Err>1

%Jacobian Evaluation
[Rvdd,Rvd,Rv,Rhs,qnm,Dnm]=JacobInd0(tnm,vnm,vdnm,vddnm,unm,par,q0,V,U,B,...
    SJDT,SMDT,STSDAT);
Phiqnm=PhiqEval(tnm,qnm,SJDT,par);
J=[Rvdd+(h/2)*Rvd+((h^2)/4)*Rv,Phiqnm'];
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

vn2=vnm+h*vdnm+((h^2)/2)*vddnm+((h^3)/6)*vdddnm;    %Solution estimate for 
vdn2=vdnm+h*vddnm+((h^2)/2)*vdddnm;                 %error control                                      

%Solution Estimates
vdd=vddnm+h*vdddnm; %Solution estimate for integration
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);
Lam=Lamnm+h*Lamdnm;
u=unm;

% Solve Discretized Equations

i=1;        %Set solution iteration counter
err=intol+1;

while err>intol

%Residual Calculation
R=ResidInd0(tn,v,vd,vdd,Lam,u,B,q0,V,U,SJDT,SMDT,STSDAT,par);

if i==1
    R1n=norm(R);
end

% Newton Correction

x=-J\R;

vdd=vdd+Pvdd*x;
Lam=Lam+PLam*x;

%Evaluate v and vd
vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);

err=norm(R);
i=i+1;

end
Jiter=i-1;

if hvar==1  %variable step

%Evaluate  error with sciv=Atol*(1+abs(v(i))), scivd=Atol*(1+abs(vd(i)))...
%and p=2
vdiff=v-vn2;
vddiff=vd-vdn2;
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

hopt=h*(1/Err)^(1/3);

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
if h<10^-5
h=10^-5;
Err=0.5;
    end
end

if h>hmax
    h=hmax; 
end

end

end

if hvar==2
Err=0.5;
end

end


end



