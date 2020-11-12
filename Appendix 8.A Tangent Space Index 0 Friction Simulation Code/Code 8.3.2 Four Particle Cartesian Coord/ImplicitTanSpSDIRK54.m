function[v,vd,vdd,Lam,jiter,R1n,JCond,h,err,nch]=ImplicitTanSpSDIRK54(n,...
    tn,npar,Vv,Vvd,Vvdd,LLam,Uu,Q,Qd,Qdd,q0,U,V,B,par,intol,Atol,nq,nh,...
    h,hmax,nch,hvar)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

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

vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);

if n-1==npar
vdd=vddnm;
Lam=LLam(:,n-1);
end

if n-1>npar
vdd=2*vddnm-Vvdd(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

vd=vdnm+(h/2)*(vddnm+vdd);
v=vnm+h*vdnm+((h^2)/4)*(vddnm+vdd);

%Jacobian Evaluation
ue=Uu(:,n-1);
[u,Iteru] = usolv(ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
[B,Biter]=CorrectB(q,B,U,par);
[QAsq,QAsqd,QAsLam]=QAsqqdeval(q,qd,Lam,par);
[Rvdd,Rvd,Rv,Phiq]=JacobFull(v,vd,vdd,Lam,u,par,q0,V,U,B);
J=[Rvdd+h*a*Rvd+(h^2)*aa*Rv,Phiq'-QAsLam];
JCond=cond(J);

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

[u,Iteru] = usolv(u,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par);
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
Pst=zeros(nh,1);
qd=D*vd-U*B*Pst;
QA=QAEval(q,qd,Lam,par);
S=zeros(nq,1);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);

R1=M*D*k1+Phiq'*Lam-(M*U*B*Gam+S+QA);

if i==1
    R11Norm=norm(R1);
end

% Quasi Newton Correction
z=-J\R1;
delk1=[z(1);z(2)];
delLam=[z(3);z(4);z(5);z(6);z(7);z(8);z(9);z(10)];
k1=k1+delk1;
Lam=Lam+delLam;
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

[u,Iteru]=usolv(ue,vnm,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par);
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
Pst=zeros(nh,1);
qd=D*vd-U*B*Pst;
QA=QAEval(q,qd,Lam,par);
S=zeros(nq,1);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);

R2=M*D*k2+Phiq'*Lam-(M*U*B*Gam+S+QA);

if i==1
    R21Norm=norm(R2);
end

% Quasi Newton Correction
z=-J\R2;
delk2=[z(1);z(2)];
delLam=[z(3);z(4);z(5);z(6);z(7);z(8);z(9);z(10)];
k2=k2+delk2;
Lam=Lam+delLam;
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

[u,Iteru]=usolv(ue,vnm,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par);
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
Pst=zeros(nh,1);
qd=D*vd-U*B*Pst;
QA=QAEval(q,qd,Lam,par);
S=zeros(nq,1);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);

R3=M*D*k3+Phiq'*Lam-(M*U*B*Gam+S+QA);

if i==1
    R31Norm=norm(R1);
end

% Quasi Newton Correction
z=-J\R3;
delk3=[z(1);z(2)];
delLam=[z(3);z(4);z(5);z(6);z(7);z(8);z(9);z(10)];
k3=k3+delk3;
Lam=Lam+delLam;
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

[u,Iteru]=usolv(ue,vnm,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par);
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
Pst=zeros(nh,1);
qd=D*vd-U*B*Pst;
QA=QAEval(q,qd,Lam,par);
S=zeros(nq,1);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);

R4=M*D*k4+Phiq'*Lam-(M*U*B*Gam+S+QA);

if i==1
    R41Norm=norm(R4);
end

% Quasi Newton Correction
z=-J\R4;
delk4=[z(1);z(2)];
delLam=[z(3);z(4);z(5);z(6);z(7);z(8);z(9);z(10)];
k4=k4+delk4;
Lam=Lam+delLam;
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

[u,Iteru]=usolv(ue,vnm,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=CorrectB(q,B,U,par);
Phiq=PhiqEval(q,par);
D=(eye(nq)-U*B*Phiq)*V;
Pst=zeros(nh,1);
qd=D*vd-U*B*Pst;
QA=QAEval(q,qd,Lam,par);
S=zeros(nq,1);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);

R5=M*D*k5+Phiq'*Lam-(M*U*B*Gam+S+QA);

if i==1
    R51Norm=norm(R5);
end

% Quasi Newton Correction
z=-J\R5;
delk5=[z(1);z(2)];
delLam=[z(3);z(4);z(5);z(6);z(7);z(8);z(9);z(10)];
k5=k5+delk5;
Lam=Lam+delLam;
err=norm(R5);
i=i+1;
end
j5=i-1;

R1=[R11Norm;R21Norm;R31Norm;R41Norm;R51Norm];
R1n=norm(R1,inf);
iter=[j1;j2;j3;j4;j5];
jiter=max(iter);

k=[k1,k2,k3,k4,k5];

%Evaluate vn, vdn, and vddn
v=vnm+h*vdnm+(h^2)*k*B1';
vd=vdnm+h*k*b';

[vdd,Lam,jodeiter,ECond]=ODEfunct(n,tn,v,vd,k5,Lam,Uu,...
    q0,V,U,B,par);

if hvar==1  %variable step

%Evaluate imbedded solution for error estimate vnhat and vdnhat
vnhat=vnm+h*vdnm+(h^2)*k*D1';
vdnhat=vdnm+h*k*d1';

%Evaluate  error with sc=Atol, and p=5
diffsqr=(v-vnhat)'*(v-vnhat)+(vd-vdnhat)'*(vd-vdnhat);
err=sqrt((1/(2*(nq-nh)*(Atol^2)))*diffsqr);

%Change step size

hopt=h*(1/err)^(1/6);

if hopt<h
h=hopt/2;
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