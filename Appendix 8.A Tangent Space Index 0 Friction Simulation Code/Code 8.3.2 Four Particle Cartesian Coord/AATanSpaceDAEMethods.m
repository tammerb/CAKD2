%Tangent Space RK DAE Exp; 4-Particle Friction-Smothed Sign-Cartesian

utol=10^-8;     %Tolerance in solving for u
Btol=10^-8;     %Convergence criteria in B iteration
intol=10^-6;     %Tolerance in solving discretized equations of motion
Atol=10^-6;
Maxv=0.7;       %Limit on magnitude of v
Maxjiter=10;      %Limit on number of implicit integration iterations
MaxJcond=100;    %Maximum J condition no in implicit integration
Maxjodeiter=10;  %Limit on number of explicit ODE iterations
MaxECond=100;    %Maximum E condition no in explicit integration

h0=0.0001;       %Initial h
h=h0;
hmax=0.01;
hvar=1;     %hvar=1, variable h;hvar=2, constant h

tfinal=4;

integ=1;  %integ=1, TanSp Trap; integ=2, TanSp SDIRK54
          %integ=3, Nystrom 4; integ=4, RKFN45

nq=10;
nh=8;
nv=nq-nh;
nu=nh;

g=9.8;
K1=2;   %Spring Constants
K2=2;
K3=2;
K4=10;
m1=2;
m2=2;
m3=6;
m4=6;
el1=5;   %Length of distance constraint
el2=7;

FR=2;   %FR=1,Rectangular Joint;FR=2,Cylindrical Joint

%Friction Data
vt=20*h0;
mud=0.16;  %Coefficient of Dynamic Friction
mus=0.2;

%Enter all parameters to be used
par=[nq;nh;nv;nu;g;m1;m2;m3;m4;K1;K2;K3;K4;el1;el2;...
    mud;mus;vt;utol;Btol;intol;h0;FR];

% Data Storage Arrays
Qd=zeros(nq,10);
Q=zeros(nq,10);
Qdd=zeros(nq,10);
Vv=zeros(nv,10);
Vvd=zeros(nv,10);
Vvdd=zeros(nv,10);
Uu=zeros(nu,10);
Uud=zeros(nu,10);
Uudd=zeros(nu,10);
LLam=zeros(nh,10);

% Initial Conditions
t(1)=0;
q0=[4;0;0;0;3;0;0;0;6.32;5];
qd0=[0;0;0;0;0;0;0;0;0;-1];

%Initial qdd and Lam

Lam=zeros(8,1);
q=q0;
qd=qd0;
M=MEval(q,par);
Phiq=PhiqEval(q,par);
Gam=GamEval(q,qd,par);
err=1;
i=1;
while err >0.00001

QA=QAEval(q,qd,Lam,par);
AccCoef=[M,Phiq';Phiq,zeros(nh,nh)];
AccRHS=[QA;-Gam];

z=AccCoef\AccRHS;

AccCoefcond(i)=cond(AccCoef);

qdd=[z(1);z(2);z(3);z(4);z(5);z(6);z(7);z(8);z(9);z(10)];
Lam=[z(11);z(12);z(13);z(14);z(15);z(16);z(17);z(18)];

err=norm(M*qdd+Phiq'*Lam-QA);
errrpt(i)=err;
i=i+1;
initirpt=i;
end

Q(:,1)=q;
Qd(:,1)=qd;
Qdd(:,1)=qdd;
LLam(:,1)=Lam;
qnorm(1)=norm(q);
qdnorm(1)=norm(qd);

RMcond(1)=0;

%Start Integration Process
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
R1nrpt(1)=0;
SolIter(1)=Maxjiter+1; 
jRepar=0;       %Counter for Reparameterization
Cr=2;
jiterrpt(1)=0;
vnormrpt(1)=Maxv+1;
nch=1;
STKWork=0;
JCondrpt(1)=MaxJcond+1;
jodeiterrpt(1)=Maxjodeiter+1;
ECondrpt(1)=MaxECond+1;

while t(n)<tfinal

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

if integ<3     %Implicit Integration

if vnormrpt(n-1)>Maxv
Cr=2;
end

if jiterrpt(n-1)>Maxjiter
Cr=2;
end

if JCondrpt(n-1)>MaxJcond
Cr=2;
end 

if Cr>1

npar=n-1;    
Crrpt(n)=Cr;
    
% Parameterization

[v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar);

u=zeros(nu,1);
Uu(:,n-1)=u;
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;

jReparrpt(n)=jRepar;

Crr=1;
end

Cr=0;

if integ==1     %Implicit Trap
    
[v,vd,vdd,Lam,jiter,R1n,h,JCond,Rnorm]=ImplicitTanSpTrap(n,npar,...
    Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,par);

EqErr(n)=Rnorm;

end      %End Trap Integ

if integ==2     %Implicit SDIRK54
    
[v,vd,vdd,Lam,jiter,R1n,JCond,h,err,nch]=ImplicitTanSpSDIRK54(n,tn,...
    npar,Vv,Vvd,Vvdd,LLam,Uu,Q,Qd,Qdd,q0,U,V,B,par,intol,Atol,nq,nh,...
    h,hmax,nch,hvar);
 
hrpt(n)=h;
end    %end Implicit SDIRK54

R1Normrpt(n)=R1n;
jiterrpt(n)=jiter;
JCondrpt(n)=JCond;
    
end      %end implicit integration 

    
if integ>2      %Start Explicit integration
    
if vnormrpt(n-1)>Maxv
Cr=2;
end

if jodeiterrpt(n-1)>Maxjodeiter
Cr=2;
end

if ECondrpt(n-1)>MaxECond
Cr=2;
end

if Cr>1    %Reparameterize

npar=n-1;    
Crrpt(n)=Cr;
    
% Parameterization

[v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar);

u=zeros(nu,1);
Uu(:,n-1)=u;
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;

jReparrpt(n)=jRepar;

end   %end reparameterization
    
if integ==3     %Explicit Nystrom4
    
[v,vd,vdd,Lam,ECond,jodeiter]=ExplicitNystrom4(n,tnm,...
    Vv,Vvd,Vvdd,LLam,Uu,U,V,B,q0,par); 
    
end       %End Explicit Nystrom4 integ

if integ==4     %ExplicitRKFN45
 
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
Lamnm=LLam(:,n-1);
[v,vd,vdd,Lam,ECond,h,nch,jodeiter]=ExplicitRKFN45(n,tnm,...
    vnm,vdnm,vddnm,Lamnm,Uu,U,V,B,q0,h,hmax,par,Atol,nq,nh,nch,hvar);
hrpt(n)=h;

end   %end RKFN45

ECondrpt(n)=ECond;
jodeiterrpt(n)=jodeiter;

end       %end explicit integration

%Process Results
Vv(:,n)=v;
Vvd(:,n)=vd;
Vvdd(:,n)=vdd;
LLam(:,n)=Lam;
vnormrpt(n)=norm(v);
vdnormrpt(n)=norm(vd);
vddnormrpt(n)=norm(vdd);
Lamnormrpt(n)=norm(Lam);

%Evaluate q

ue=Uu(:,n-1);
[u,Iteru]=usolv(u,v,q0,V,U,B,par);
Iterurpt(n)=Iteru;
Uu(:,n)=u;
unorm(n)=norm(u);
q=q0+V*v-U*u;
Q(:,n)=q;
qnormrpt(n)=norm(q);

%Update B
[B,Biter]=CorrectB(q,B,U,par);
Biterrpt(n)=Biter;

% Evaluate qd
Phiq = PhiqEval(q,par);

ondPhiq(n)=cond(Phiq);
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

% Evaluate qdd

[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

%Calculate Total Energy
M=MEval(q,par);
TE(n)=0.5*qd'*M*qd+m3*g*q(9)+0.5*K4*(q(10)-q(1)-1)^2+...
    0.5*K1*q(1)^2+0.5*K2*q(5)^2+0.5*K3*q(9)^2;

%Data of Interest
x1(n)=q(1);
y2(n)=q(5);
z3(n)=q(9);
x4(n)=q(10);
x1d(n)=qd(1);
y2d(n)=qd(5);
z3d(n)=qd(9);
x4d(n)=qd(10);
x1dd(n)=qdd(1);
y2dd(n)=qdd(5);
z3dd(n)=qdd(9);
x4dd(n)=qdd(10);

%Friction Forces

fr1(n)=-mud*abs(Lam(1))*csign(qd(2),par);
fr2(n)=-mud*abs(Lam(2))*csign(qd(3),par);

if max(abs(qd(2)),abs(qd(3)))<vt
STKWork=STKWork+abs(qd(2)*fr1(n))*h+abs(qd(3)*fr2(n))*h;
STKWorkRpt(n)=STKWork;
end

qd23norm(n)=norm([qd(2);qd(3)]);
qdd23norm(n)=norm([qdd(2);qdd(3)]);

%Calculate constraint error
PosConstrNorm(n) = norm(PhiEval(q,par));
VelConstrNorm(n)=norm(Phiq*qd);
AccelErr=Phiq*qdd+Gam;
AccConstrNorm(n)=norm(AccelErr);

end


   