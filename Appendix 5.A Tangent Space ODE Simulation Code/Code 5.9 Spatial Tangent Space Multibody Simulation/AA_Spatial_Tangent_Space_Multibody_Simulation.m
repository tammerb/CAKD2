%AA Spatial Tangent Space Multibody Simulation

utol=10^-8;     %Tolerance in solving for u
Btol=10^-8;     %Convergence criteria in B iteration
intol=10^-6;     %Tolerance in solving discretized equations of motion
Atol=10^-5;     %Absolute error tolerance for variable step metnods
Maxv=10;       %Limit on magnitude of v
MaxIntiter=12;      %Limit on number of integration iterations
MaxUiter=8;      %Limit on number of U iterations
MaxJcond=200;   %Limit on magnitude of Econd in Trap
R1nmax=300;     %Limit on residual in integrationt
MaxECond=50;  %Limit on magnitude of ECond
MaxBCond=6000;  %Limit on magnitude of BCond
MaxBnormRat=10;  %Maximum ratio of Bnorm to that after parameterizationm
MaxBCondRat=10;     %Maximum ratio of BCond to BCond0 at parameterization

h0=0.001;       %Initial time step
h=h0;
hmax=0.001;
hvar=1;         %hvar=1, variable h;hvar=2, constant h

tfinal=2;

integ=1;            %integ=1-Impl ODE Trap; integ=2-Impl ODE SDIRK54; 
                    %integ=3-Impl Ind0 Trap; integ=4-Impl Ind0 SDIRK54;
                    %integ=5-Expl ODE Nystrom4; integ=6-Expl ODE RKFN45;
                    %integ=7-Expl Ind0 Nystrom4; integ=8-Expl Ind0 RKFN45
            
InvJ=2;      %In Implicit Integration, InvJ=1, Invert J; InvJ=2, Factor J            
            
g=9.8;

%Application Data
    %app=1, Spin Stabilized Top
    %app=2, Spatial Double Pendulum
    %app=3, One Body Cylindrical with Spring
    %app=4, Spatial Slider-Crank
    %app=5, 4-Translating Mass Model
    %app=6, Rotating Disk with Translating Body
    
app=2;

[nb,ngc,nh,nc,nv,nu,NTSDA,SJDT,SMDT,STSDAT,q0,qd0]=AppData(app);

par=[nb;ngc;nh;nc;nv;nu;g;utol;Btol;intol;Atol;h0;hvar;NTSDA];

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

Q(:,1)=q0;
Qd(:,1)=qd0;

%Calculate Initial Acceleration and Lagrange Multipliers
M=MEval(q0,SMDT,par);
G=GEval(q0);
Phiq0=PhiqEval(0,q0,SJDT,par);

Phiq0Cond=cond(Phiq0);

QA=QAEval(0,q0,qd0,SMDT,STSDAT,par);
S=SEval(q0,qd0,SMDT,par);
Gam=GamEval(0,q0,qd0,SJDT,par);
EE=[M,Phiq0';Phiq0,zeros(nc,nc)];
EEcond=cond(EE);
RHS=[QA+S;-Gam];

x=EE\RHS;

Pqdd=[eye(ngc),zeros(ngc,nc)];
PLam=[zeros(nc,ngc),eye(nc)];
qdd=Pqdd*x;
Lam=PLam*x;

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

if integ<5     %Start implicit reparameterization criteria
    
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

if integ>4      %Start explicit reparameterization criteria
    
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

if integ==1     %Implicit ODE Trap    
[v,vd,vdd,R1n,Jiter,JCond,Jinv,Jinviter,h,nch]=...
    ImplicitODETrap(n,tn,npar,...
    Vv,Vvd,Vvdd,Uu,q0,V,U,B,h,hmax,nch,SMDT,STSDAT,SJDT,par,InvJ,Jinv);

R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
Jiterrpt(n)=Jiter;
hrpt(n)=h;
end

if integ==2     %Implicit ODE SDIRK54
[v,vd,vdd,Jiter,R1Norm,J,JCond,h,nch,Jinv,Jinviter]=...
    ImplicitODESDIRK54(n,tn,Vv,Vvd,Vvdd,Uu,q0,U,V,B,par,...
    h,hmax,nch,npar,SJDT,SMDT,STSDAT,Jinv,InvJ);

R1Normrpt(n)=R1Norm;
JCondrpt(n)=JCond;
Jiterrpt(n)=Jiter;
hrpt(n)=h;
end

if integ==3     %Implicit Ind0 Trap   
[v,vd,vdd,Lam,R1n,Jiter,JCond,h,nch]=...
    ImplicitInd0Trap(n,tn,npar,Vv,Vvd,Vvdd,LLam,Uu,q0,V,U,B,h,hmax,nch,...
    SMDT,STSDAT,SJDT,par); 

R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
Jiterrpt(n)=Jiter;
hrpt(n)=h;
LLam(:,n)=Lam;
end

if integ==4     %Implicit Ind0 SDIRK54 
[v,vd,vdd,Lam,Jiter,R1n,J,JCond,h,nch]=...
    ImplicitInd0SDIRK54(n,tn,Vv,Vvd,Vvdd,LLam,Uu,q0,U,V,B,par,...
    h,hmax,nch,npar,SJDT,SMDT,STSDAT);

R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
Jiterrpt(n)=Jiter;
hrpt(n)=h;
LLam(:,n)=Lam;
end

if integ==5     %Explicit ODE Nystrom4    
[v,vd,vdd,ECond]=ExplicitODENystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,SMDT,STSDAT,SJDT,par);
ECondrpt(n)=ECond;    
end

if integ==6     %Explicit ODE RKFN45
[v,vd,vdd,ECond,h,nch]=ExplicitODERKFN45(n,tn,Vv,Vvd,Uu,...
    U,V,B,q0,h,hmax,par,SMDT,STSDAT,SJDT,nch);

ECondrpt(n)=ECond;
hrpt(n)=h;
end

if integ==7     %Explicit Ind0 Nystrom4
[v,vd,vdd,Lam,ECond]=ExplicitInd0Nystrom4(n,tn,Vv,Vvd,Uu,...
    V,U,B,q0,h,SMDT,STSDAT,SJDT,par);

ECondrpt(n)=ECond;
LLam(:,n)=Lam;
end

if integ==8     %Explicit Ind0 RKFN45
[v,vd,vdd,Lam,ECond,h,nch]=ExplicitInd0RKFN45(n,tn,Vv,Vvd,Uu,...
    U,V,B,q0,h,hmax,par,SMDT,STSDAT,SJDT,nch);

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
BCond(n)=Bnorm(n)/norm(Phiq*U);
BCondRat(n)=BCond(n)/BCond0;
Binvnorm(n)=norm(Phiq*U);
Bcond(n)=Bnorm(n)*Binvnorm(n);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

% Evaluate qdd and Lam (Lam evaluation a postprocessing step)

Gam=GamEval(tn,q,qd,SJDT,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

Lam=B'*U'*(-MEval(q,SMDT,par)*qdd+QAEval(tn,q,qd,SMDT,STSDAT,par));
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
T=1;
while T<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=STSDATPart(STSDAT,T);
[r1,p1]=qPart(q,i);
r2=[0;0;0];
p2=[1;0;0;0];
r1d=[0;0;0];
p1d=zeros(4,1);
if j>=1
[r2,p2]=qPart(q,j);
end
A1=ATran(p1);
A2=ATran(p2);
d12=r2+A2*s2pr-r1-A1*s1pr;
el=sqrt(d12'*d12);
SET=0.5*K*(el-el0)^2;
SE(n)=SE(n)+SET;
T=T+1;
end

TE(n)=KE(n)+PE(n)+SE(n);

%Data of Interest (Enter for each application)

if app==1   %Spin Stabilized Top
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Spatial Double Pendulum
x2(n)=q(8);
y2(n)=q(9);
z2(n)=q(10);
x2d(n)=qd(8);
y2d(n)=qd(9);
z2d(n)=qd(10);
x2dd(n)=qdd(8);
y2dd(n)=qdd(9);
z2dd(n)=qdd(10);
Lam4(n)=Lam(4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Body in Cylindrical Joint with Spring
[r,p]=qPart(q,1);
[rd,pd]=qPart(qd,1);
z(n)=r'*uz;
zd(n)=rd'*uz;
omegz(n)=2*uz'*EEval(p)*pd;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Spatial Slider-Crank
x2(n)=q(8);
x2d(n)=qd(8);
x2dd(n)=qdd(8);
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
omeg1x(n)=(2*EEval(p1)*p1d)'*ux;
omeg1xd(n)=(2*EEval(p1)*p1dd)'*ux;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5   %Four Translating Mass Spatial Model
QA=QAEval(tn,q,qd,SMDT,STSDAT,par);
QARec(:,n)=QA;
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
if app==6       %Rotating Disk with Translating Body
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate constraint error
Phi=PhiEval(tn,q,SJDT,par);
Phiq=PhiqEval(tn,q,SJDT,par);
%[Gam,Gamsq,Gamsqd]=GamEval(tn,q,qd,SJDT,par);
PosConstrNorm(n)=norm(Phi);
VelConstrNorm(n)=norm(Phiq*qd);
AccConstrNorm(n)=norm(Phiq*qdd+Gam);

end



