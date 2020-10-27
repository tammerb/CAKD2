%AA Spatial Index 0 Multibody Simulation with Friction

utol=10^-8;     %Tolerance in solving for u
Btol=10^-8;     %Convergence criteria in B iteration
intol=10^-6;     %Tolerance in solving discretized equations of motion
Atol=10^-5;     %Absolute error tolerance for variable step metnods
Vtol=10^-4;         %Velocity error tolerance for variable step metnods
Maxv=10;       %Limit on magnitude of v
MaxIntiter=12;      %Limit on number of integration iterations
MaxUiter=8;      %Limit on number of U iterations
MaxJcond=100;   %Limit on magnitude of Econd in Trap
R1nmax=30000;     %Limit on residual in integrationt
MaxECond=20;  %Limit on magnitude of ECond
MaxBCond=6000;  %Limit on magnitude of BCond
MaxBnormRat=10;  %Maximum ratio of Bnorm to that after parameterizationm
MaxBCondRat=10;     %Maximum ratio of BCond to BCond0 at parameterization

h0=0.0001;       %Initial time step
h=h0;
hmax=0.001;
hvar=1;     %hvar=1, variable h;hvar=2, constant h

tfinal=0.25;

integ=3;            %integ=1-Impl Ind0 Trap; integ=2-Impl Ind0 SDIRK54; 
                    %integ=3-Expl Ind0 Nystrom4; integ=4-Expl Ind0 RKFN45          
            
g=9.8;              %gravitational acceleration
vt=100*h0;           %Transition time in continuous friction function

%Application Data
    %app=1, One Body Cylindrical Joint with Ground and Spring
    %app=2, Spatial Slider-Crank
    %app=3, Four Body Slider
    %app=4, %Rotating Disk with Translating Body
    %app=5, %One Body Pendulum, Cyl about x axis
    %app=6, %One Body Translation Along x axis
    
app=2;

[nb,ngc,nh,nc,nv,nu,NTSDA,SJDT,SMDT,STSDAT,q0,qd0]=AppData(app);

par=[nb;ngc;nh;nc;nv;nu;g;utol;Btol;intol;Atol;Vtol;hvar;NTSDA;vt];

if 1>2
%Check initial conditions for feasibility
Phi0=PhiEval(0,q0,SJDT,par);
if norm(Phi0)>10^-3
infeasableconfig=1;
%return
end
Phiq0=PhiqEval(0,q0,SJDT,par);
Phid0=Phiq0*qd0;
if norm(Phid0)>10^-3
infeasablevelocity=1; 
%return
end
end

%Initial condition calculation, if required

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

% Data Storage Arrays
Q=zeros(ngc,10);
Qd=zeros(ngc,10);
Qdd=zeros(ngc,10);
Vv=zeros(nv,10);
Vvd=zeros(nv,10);
Vvdd=zeros(nv,10);
Uu=zeros(nu,10);
Uud=zeros(nu,10);
Uudd=zeros(nu,10);
LLam=zeros(nc,10);
LLam0=zeros(nc,10);
QQA=zeros(ngc,10);

Q(:,1)=q0;
Qd(:,1)=qd0;

%Calculation of initial values of qdd and Lam to start integration
%Coefficients of friction are indexed by w/N, w=0,1,...N in Function 
%Ind0IC to calculate a sequence of qdd and Lam that converge to the
%desired coefficients of friction. The user may select N.

N=10;
[qdd,Lam,Qdd0,LLam0,w]=Ind0IC(q0,qd0,SMDT,SJDT,STSDAT,par,N);

Qdd(:,1)=qdd;
LLam(:,1)=Lam;

%Initialize Data For Integration
    
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
jRepar=0;       %Counter for Reparameterization
Cr=2;
Jiterrpt(1)=0;
Jiter=0;
Uterurpt(1)=0;
JCondrpt(1)=0;
R1Normrpt(1)=R1nmax+1;
vnormrpt(1)=Maxv+1;
ECondrpt(1)=MaxECond+1;
Econd(1)=1;
uiterrpt(1)=1;
nch=1;
VelConstrNorm(1)=0;
AccConstrNorm(1)=0;
nres=1;
jnres=0;
B0norm=1;
BnormRat(1)=1;
BCond(1)=MaxBCond+1;
BCondRat(1)=1;
Jinv=eye(nv);

%Integration
while t(n)<tfinal

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

if integ<3     %Start implicit reparameterization criteria
    
if vnormrpt(n-1)>Maxv
Cr=2;
end

if Jiterrpt(n-1)>MaxIntiter
Cr=Cr+2;
end

if JCondrpt(n-1)>MaxJcond
Cr=Cr+2;
end

if R1Normrpt(n-1)>R1nmax
    Cr=Cr+2;
end

if BnormRat(n-1)>MaxBnormRat
    Cr=Cr+2;
end

end      %end implicit reparameterization criteria

if integ>2      %Start explicit reparameterization criteria
    
if vnormrpt(n-1)>Maxv
Cr=2;
end
    
if ECondrpt(n-1)>MaxECond
Cr=Cr+2;
end

if BnormRat(n-1)>MaxBnormRat
    Cr=Cr+2;
end
    
end            %end explicit reparameterization criteria

if Cr>1

npar=n-1;    
Crrpt(n)=Cr;
nch=1;
    
% Parameterization

tnm=t(n-1);

[v,vd,vdd,q0,U,V,B,jRepar]=Param(n,tnm,Q,Qd,Qdd,SJDT,par,jRepar);

u=zeros(nu,1);
Uu(:,n-1)=u;
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;

B0norm=norm(B);
BCond0=B0norm*norm(U'*U);

jReparrpt(n)=jRepar;

Cr=1;
end

% Integration

if integ==1     %Implicit Ind0 Trap   
[v,vd,vdd,Lam,R1n,Jiter,JCond,h,nch]=...
    ImplicitInd0Trap(n,tn,npar,Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,h,hmax,nch,...
    SMDT,SJDT,STSDAT,par);

R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
Jiterrpt(n)=Jiter;
hrpt(n)=h;
LLam(:,n)=Lam;
end

if integ==2     %Implicit Ind0 SDIRK54 
[v,vd,vdd,Lam,Jiter,R1n,J,JCond,h,nch]=...
    ImplicitInd0SDIRK54(n,tn,Vv,Vvd,Vvdd,LLam,Uu,q0,U,V,B,par,...
    h,hmax,nch,npar,SJDT,SMDT,STSDAT);

R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
Jiterrpt(n)=Jiter;
hrpt(n)=h;
LLam(:,n)=Lam;
end

if integ==3     %Explicit Ind0 Nystrom4
[v,vd,vdd,Lam,ECond,jodeiter1]=ExplicitInd0Nystrom4(n,tn,Vv,Vvd,...
    Vvdd,LLam,Uu,V,U,B,q0,h,npar,SJDT,SMDT,STSDAT,par);

ECondrpt(n)=ECond;
LLam(:,n)=Lam;
jodeiterRpt(n)=jodeiter1;
end

if integ==4     %Explicit Ind0 RKFN45
[v,vd,vdd,Lam,ECond,h,nch]=ExplicitInd0RKFN45(n,tn,Vv,Vvd,...
    Vvdd,LLam,Uu,U,V,B,q0,h,hmax,par,SMDT,STSDAT,SJDT,nch);

ECondrpt(n)=ECond;
LLam(:,n)=Lam;
hrpt(n)=h;
end

%Process Results
Vv(:,n)=v;
Vvd(:,n)=vd;
Vvdd(:,n)=vdd;
vnormrpt(n)=norm(v);
vdnormrpt(n)=norm(vd);
vddnormrpt(n)=norm(vdd);

%Evaluate q
u=Uu(:,n-1);
[u,Iteru]=usolv(tn,u,v,q0,SJDT,V,U,B,par);
Iterurpt(n)=Iteru;
Uu(:,n)=u;
unorm(n)=norm(u);
q=q0+V*v-U*u;
Q(:,n)=q;
qnormrpt(n)=norm(q);

%Update B and Evaluate qd
[B,Biter]=CorrectB(tn,q,B,U,SJDT,par);
Biterrpt(n)=Biter;
Bnorm(n)=norm(B);
BnormRat(n)=Bnorm(n)/B0norm;
Phiq=PhiqEval(tn,q,SJDT,par);
BCond(n)=Bnorm(n)*norm(Phiq*U);
BCondRat(n)=BCond(n)/BCond0;
Binvnorm(n)=norm(Phiq*U);
Bcond(n)=Bnorm(n)*Binvnorm(n);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

% Evaluate qdd

Gam=GamEval(tn,q,qd,SJDT,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);
LLam(:,n)=Lam;

%End Tangent Space Integration

%Calculate Total Energy
M=MEval(q,SMDT,par);
KE(n)=0.5*qd'*M*qd;

i=1;
PE(n)=0;
while i<=nb
 PE(n)=PE(n)+SMDT(1,i)*g*q(7*(i-1)+3);
 i=i+1;
end

SE(n)=0;
TS=1;
while TS<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=STSDATPart(STSDAT,TS);
[r1,p1]=qPart(q,i);
r2=[0;0;0];
p2=[1;0;0;0];
if j>=1
[r2,p2]=qPart(q,j);
end
A1=ATran(p1);
A2=ATran(p2);
d12=r2+A2*s2pr-r1-A1*s1pr;
el=sqrt(d12'*d12);
SET=0.5*K*(el-el0)^2;
SE(n)=SE(n)+SET;
TS=TS+1;
end

TE(n)=KE(n)+PE(n)+SE(n);

%Data of Interest (Enter for each application)

if app==1   %Body in Cylindrical Joint with Ground
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
z1d(n)=qd(3);
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
omegaz(n)=uz'*2*EEval(p1)*p1d;
[fx1prk,fy1prk,fx2prk,fy2prk,f12kcylfr,tau12kcylfr]=...
    QAcylfEval(tn,q,qd,Lam,SMDT,SJDT,STSDAT,par);
fx1prkrpt(n)=fx1prk;
fy1prkrpt(n)=fy1prk;
fx2prkrpt(n)=fx2prk;
fy2prkrpt(n)=fy2prk;
f12kcylfrrpt(n)=f12kcylfr;
tau12kcylfrrpt(n)=tau12kcylfr;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Spatial Slider-Crank
x2(n)=q(8);
x2d(n)=qd(8);
x2dd(n)=qdd(8);
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
omeg1x(n)=(2*EEval(p1)*p1d)'*ux;
omeg1xd(n)=(2*EEval(p1)*p1dd)'*ux;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Four Body Slider
x1(n)=q(1);
y2(n)=q(9);
z3(n)=q(17);
x4(n)=q(22);
x1d(n)=qd(1);
y2d(n)=qd(9);
z3d(n)=qd(17);
x4d(n)=qd(22);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4       %Rotating Disk with Translating Body
theta1(1)=0;
p1=[q(4);q(5);q(6);q(7)];
p1d=[qd(4);qd(5);qd(6);qd(7)];
r2=[q(8);q(9);q(10)];
E1=EEval(p1);
A1=ATran(p1);
omegaz1(n)=2*uz'*E1*p1d;
theta1(n)=theta1(n-1)+h*omegaz1(n);
dely2pr(n)=(A1'*(r2-A1*ux))'*uy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5   %One Body Pendulum, Cyl
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
omeg1x(n)=ux'*2*EEval(p1)*p1d;
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
x1d(n)=qd(1);
y1d(n)=qd(2);
z1d(n)=qd(3);
x1dd(n)=qdd(1);
y1dd(n)=qdd(2);
z1dd(n)=qdd(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==6   %One Body Translation Along x Axis
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
omeg1x(n)=ux'*2*EEval(p1)*p1d;
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
x1d(n)=qd(1);
y1d(n)=qd(2);
z1d(n)=qd(3);
x1dd(n)=qdd(1);
y1dd(n)=qdd(2);
z1dd(n)=qdd(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate constraint error
Phi=PhiEval(tn,q,SJDT,par);
Phiq=PhiqEval(tn,q,SJDT,par);
%[Gam,Gamsq,Gamsqd]=GamEval(tn,q,qd,SJDT,par);
PosConstrNorm(n)=norm(Phi);
VelConstrNorm(n)=norm(Phiq*qd);
AccConstrNorm(n)=norm(Phiq*qdd+Gam);

end



