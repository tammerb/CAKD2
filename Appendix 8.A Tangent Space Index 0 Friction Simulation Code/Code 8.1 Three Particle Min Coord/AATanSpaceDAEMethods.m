%Tangent Space RK DAE Exp; 3-Particle Friction-Smothed Sign

utol=10^-10;     %Tolerance in solving for u
Btol=10^-10;     %Convergence criteria in B iteration
intol=10^-8;     %Tolerance in solving discretized equations of motion
Maxv=0.7;       %Limit on magnitude of v
Maxiter=6;      %Limit on number of integration iterations

h0=0.0001;      %Initial value of h
h=h0;
tfinal=4;

integ=2;    %integ=1, TanSp Trap; integ=2, TanSp Nystrom4

nq=3;
nh=1;
nv=nq-nh;
nu=nh;

g=9.8;
K1=10;   %Spring Constants
K2=10;
m1=6.3461;
m2=2;
m3=2;

el=5;   %Length of distance constraint
mud=0.3;  %Coefficient of Dynamic Friction
mus=0.5;
vt=20*h0;

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
t0=t(1);
q0=[0;5;6];
qd0=[0;0;-1];

%Initial qdd and Lam

    Lam0=0;

Phiq=PhiqEval(q0,par);
Gam=GamEval(q0,qd0,par);

QA=QAEval(t0,q0,qd0,Lam0,par);
AccCoef=[MEval(q0,par),Phiq';Phiq,0];
AccRHS=[QA;-Gam];

z=AccCoef\AccRHS;

AccCoefcond=cond(AccCoef);

qdd0=[z(1);z(2);z(3)];
Lam0=z(4);

Q(:,1)=q0;
Qd(:,1)=qd0;
Qdd(:,1)=qdd0;
LLam(:,1)=Lam0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);

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

if vnormrpt(n-1)>Maxv
Cr=2;
end

if integ ==1    
jiter=jiterrpt(n-1);   
if jiter>Maxiter
Cr=2;
end
end

if Cr>1

npar=n-1;    
Crrpt(n)=Cr;
    
% Parameterization

[v,vd,vdd,Lam,q0,U,V,B,jRepar]=Param(n,tn,Q,Qd,Qdd,LLam,Uu,...
    par,jRepar);

u=zeros(nu,1);
Uu(:,n-1)=u;
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;

jReparrpt(n)=jRepar;

end

Cr=0;

% Integration

if integ ==1    % Tan Sp Trap Integration
[v,vd,vdd,Lam,jiter,R1n,h,JCond,Rnorm]=IntegTanSpTrap(tn,...
   n,npar,Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,par);

jiterrpt(n)=jiter;
EqErr(n)=Rnorm;
R1nrpt(n)=R1n;
JCondrpt(n)=JCond;

end           %End Tan Sp Trap Integration

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

% Evaluate qdd

[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

%Calculate Total Energy
M=MEval(q,par);
TE(n)=0.5*qd'*M*qd+m1*g*q(1)+0.5*K2*(q(3)-q(2)-1)^2+0.5*K1*q(1)^2;

%Data of Interest
y(n)=q(1);
x1(n)=q(2);
x2(n)=q(3);
yd(n)=qd(1);
x1d(n)=qd(2);
x2d(n)=qd(3);
ydd(n)=qdd(1);
x1dd(n)=qdd(2);
x2dd(n)=qdd(3);

qd12norm(n)=norm([qd(1);qd(2)]);
qdd12norm(n)=norm([qdd(1);qdd(2)]);

%Calculate constraint error
Phiq=PhiqEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd,par);
PosConstrNorm(n) = norm(PhiEval(q,par));
VelConstrNorm(n)=norm(Phiq*qd);
AccelErr=Phiq*qdd+Gam;
AccConstrNorm(n)=norm(AccelErr);

end


   