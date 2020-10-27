function[vn,vdn,vddn,Maxjiter,R1Norm,JCond,h,err]=ImplicitSDIRK54(n,tn,...
    Vv,Vvd,Vvdd,Uu,Q,Qd,Qdd,q0,U,V,B,par,intol,Atol,nq,nh,h,hmax,nch,hvar)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

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

Inv=1;    % Inv=1, factor each use; Inv=2, Use Inverse of J

% Integration Jacobian Evaluation
q=Q(:,n-1);
qd=Qd(:,n-1);
qdd=Qdd(:,n-1);
[Rv,Rvd,Rvdd]=Jacob(tn,q,qd,qdd,U,V,B,par);
J=Rvdd+h*a*Rvd+(h^2)*aa*Rv;
JCond=cond(J);
if Inv==2
    Jinv=inv(J);
end
I4=eye(4);

% Start Integration Step
q=Q(:,n-1);
qd=Qd(:,n-1);
qdd=Qdd(:,n-1);
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
ue=Uu(:,n-1);

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

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);

R1=D'*M*D*k1-D'*(M*U*B*Gamma+S+QA);

if i==1
    R11Norm=norm(R1);
end

% Quasi Newton Correction
if Inv==1
delk1=-J\R1;
end
if Inv==2
delk1=-Jinv*R1;
end
k1=k1+delk1;
err=norm(R1);
i=i+1;
end
j1=i-1;
ue=u;

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

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);

R2=D'*M*D*k2-D'*(M*U*B*Gamma+S+QA);

if i==1
    R21Norm=norm(R2);
end

% Quasi Newton Correction
if Inv==1
delk2=-J\R2;
end
if Inv==2
delk2=-Jinv*R2;
end

k2=k2+delk2;
err=norm(R2);
i=i+1;
end
j2=i-1;
ue=u;


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

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);

R3=D'*M*D*k3-D'*(M*U*B*Gamma+S+QA);

if i==1
    R31Norm=norm(R1);
end

% Quasi Newton Correction
if Inv==1
delk3=-J\R3;
end
if Inv==2
delk3=-Jinv*R3;
k3=k3+delk3;
err=norm(R3);
i=i+1;
end
j3=i-1;
ue=u;


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

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);

R4=D'*M*D*k4-D'*(M*U*B*Gamma+S+QA);

if i==1
    R41Norm=norm(R4);
end

% Quasi Newton Correction
if Inv==1
delk4=-J\R4;
end
if Inv==2
delk4=-Jinv*R4;
end
k4=k4+delk4;
err=norm(R4);
i=i+1;
end
j4=i-1;
ue=u;


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

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);

R5=D'*M*D*k5-D'*(M*U*B*Gamma+S+QA);

if i==1
    R51Norm=norm(R5);
end

% Quasi Newton Correction
if Inv==1
delk5=-J\R5;
end
if Inv==2
delk5=-Jinv*R5;
end
k5=k5+delk5;
err=norm(R5);
i=i+1;
end
j5=i-1;

R1=[R11Norm;R21Norm;R31Norm;R41Norm;R51Norm];
R1Norm=norm(R1,inf);
jiter=[j1;j2;j3;j4;j5];
Maxjiter=norm(jiter,inf);

k=[k1,k2,k3,k4,k5];

%Evaluate vn, vdn, and vddn
vn=vnm+h*vdnm+(h^2)*k*B1';
vdn=vdnm+h*k*b';

[f,E,B]=ODEfunct(tn,vn,vdn,ue,U,V,B,q0,par);
vddn=f;

if hvar==1  %variable step

%Evaluate imbedded solution for error estimate vnhat and vdnhat
vnhat=vnm+h*vdnm+(h^2)*k*D1';
vdnhat=vdnm+h*k*d1';

%Evaluate  error with sc=Atol, and p=5
diffsqr=(vn-vnhat)'*(vn-vnhat)+(vdn-vdnhat)'*(vdn-vdnhat);
err=sqrt((1/(2*(nq-nh)*(Atol^2)))*diffsqr);

%Change step size

hopt=h*(1/err)^(1/6);

if hopt<h
h=0.5*hopt;
nch=n;
end
if hopt>h
    if n>nch+5
h=hopt;
nch=n;
    end
end

if h>hmax
    h=hmax; 
end

end

end