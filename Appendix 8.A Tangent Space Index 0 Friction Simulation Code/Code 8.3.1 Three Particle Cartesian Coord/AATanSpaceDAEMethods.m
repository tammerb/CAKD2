%Tangent Space RK DAE Exp; 3-Particle Friction-Smothed-Cartesian

utol=10^-8;     %Tolerance in solving for u
Btol=10^-8;     %Convergence criteria in B iteration
intol=10^-6;     %Tolerance in solving discretized equations of motion
Maxv=0.7;       %Limit on magnitude of v
Maxiter=6;      %Limit on number of integration iterations

h0=0.0001;      %Initial h
h=h0;
tfinal=3.99;

integ=1;        %integ=1, Trap; integ=2, Nystrom4

nq=5;
nh=3;
nv=nq-nh;
nu=nh;

g=9.8;
K1=10;   %Spring Constants
K2=10;
m1=5;
m2=2;
m3=2;
el=5;   %Length of distance constraint

mud=0.08;    %Coefficient of Dynamic Friction
mus=0.1;    %Coefficient of Dynamic Friction
vt=20*h0;   %Half rise time

%Enter all parameters to be used
par=[nq;nh;nv;nu;g;m1;m2;m3;K1;K2;el;mud;mus;vt;utol;Btol;intol;h0];

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
q0=[0;0;5;0;6];
qd0=[0;0;0;0;-1];

%Initial qdd and Lam

Lam0=zeros(3,1);

Phiq=PhiqEval(q0,par);
Gam=GamEval(q0,qd0,par);

QA=QAEval(q0,qd0,Lam0,par);

AccCoef=[MEval(q0,par),Phiq';Phiq,zeros(3,3)];
AccRHS=[QA;-Gam];

z=AccCoef\AccRHS;

AccCoefcond=cond(AccCoef);

qdd0=[z(1);z(2);z(3);z(4);z(5)];
Lam0=[z(6);z(7);z(8)];

Q(:,1)=q0;
Qd(:,1)=qd0;
Qdd(:,1)=qdd0;
LLam(:,1)=Lam0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);
q=q0;
qd=qd0;
qdd=qdd0;

RMcond(1)=0;

%Start Integration Process
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
R1nrpt(1)=0;
SolIter(1)=Maxiter+1; 
jRepar=0;       %Counter for Reparameterization
Cr=2;
jiterrpt(1)=0;
vnormrpt(1)=Maxv+1;
STKWork=0;

while t(n)<tfinal

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

jiter=jiterrpt(n-1);

if vnormrpt(n-1)>Maxv
Cr=2;

end

if jiter>Maxiter
Cr=Cr+2;
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

% Integration

if integ==1
[v,vd,vdd,Lam,jiter,R1n,h,JCond,Rnorm]=ImplicitTanSpTrap(n,npar,...
    Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,par);

EqErr(n)=Rnorm;
R1nrpt(n)=R1n;
JCondrpt(n)=JCond;
jiterrpt(n)=jiter;
end

if integ==2
[v,vd,vdd,Lam,ECond,jodeiter]=ExplicitNystrom4(n,tnm,...
    Vv,Vvd,Vvdd,LLam,Uu,U,V,B,q0,par);

ECondrpt(n)=ECond;
jodeiterrpt(n)=jodeiter;
end


%Process Results
Vv(:,n)=v;
Vvd(:,n)=vd;
Vvdd(:,n)=vdd;
LLam(:,n)=Lam;
vnormrpt(n)=norm(v);
vdnormrpt(n)=norm(vd);
vddnormrpt(n)=norm(vdd);
Lamnormrpt(n)=norm(Lam);

jiterrpt(n)=jiter;

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
D=(eye(nq)-U*B*Phiq)*V;
qd=D*vd;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

%Calculate Total Energy
M=MEval(q,par);
TE(n)=0.5*qd'*M*qd+m1*g*q(2)+0.5*K2*(q(5)-q(3)-1)^2+0.5*K1*q(2)^2;

% Evaluate qdd

[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

%Data of Interest
y1(n)=q(2);
x2(n)=q(3);
x3(n)=q(5);
y1d(n)=qd(2);
x2d(n)=qd(3);
x3d(n)=qd(5);
y1dd(n)=qdd(2);
x2dd(n)=qdd(3);
x3dd(n)=qdd(5);


%Friction Forces

fr1(n)=-mud*abs(Lam(1))*sign(qd(2));
fr2(n)=-mud*abs(Lam(2))*sign(qd(3));

N1(n)=Lam(1);
N2(n)=Lam(2);

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


   