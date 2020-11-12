%AA_Planar_Multibody_Sim_DAE

%Integration and Error Control Parameters
intol=10^-6;        %Tolerance in solving discretized equations of motion
Atol=10^-4;         %Absolute error tolerance for variable step metnods
MaxIntiter=8;       %Limit on number of implicit integration iterations
MaxJcond=100;       %Limit on magnitude of Jcond in Trap
R1nmax=1000;        %Limit on residual in integrationt
MaxECond=20;        %Limit on magnitude of ECond
PosConstrMax=10^-4;  %Limit on position constraint error
VelConstrMax=10^-3;  %Limit on velocity constraint error
AccConstrMax=10^-2;  %Limit on acceleration constraint error

h0=0.001;           %Initial time step
h=h0;
hmax=0.001;         %Maximum time step

hvar=1;             %hvar=1, variable h; hvar=2, constant h
constrcor=2;         %constrcor=1, correct; constrcor=2, no correct

tfinal=10;           %Final simulation time

integ=3;    %integ=1-ImplicitIndex1; alpha=0-Trapezoidal; alpha<0-HHt
            %integ=2-ImplicitIndex2; alpha=0-Trapezoidal; alpha<0-HHt
            %integ=3-ImplicitIndex3; alpha=0-Trapezoidal; alpha<0-HHt
            %integ=4-ExplicitRKFN45
            %integ=5-ExplicitNystrom4
            
alpha=-1/3;    %Enter -1/3<=alpha<0 for HHT            
            %alpha not defined for ExplicitNystrom4 or ExplicitRKFN45
            
g=9.8;

%Application Data
    %app=1, Quick Return
    %app=2, Loader
    %app=3, Slider-Crank
    %app=4, Multiple Slider-Crank
    %app=5, Pendulum-Distance constraint
    %app=6, Pendulum-Revolute constraint
    %app=7, Double Pendulum-Unilateral spring
    %app=8, Double Pendulum
    %app=9, Lumped Mass Coil Spring-5 masses
    %app=10, Lumped Mass Coil Spring-10 masses
    
app=7;

[nb,ngc,nh,nc,NTSDA,NRSDA,PJDT,PMDT,PTSDAT,PRSDAT,q,qd]=AppData(app);    

par=[nb;ngc;nh;nc;g;intol;Atol;h0;hvar;NTSDA;NRSDA];

%Initial condition calculation, if required
if app==1       
phi2d0=6;
Phiq=PhiqEval(0,q,PJDT,par);
EEE=[Phiq;zeros(1,5),1,zeros(1,6)];
RRHS=[zeros(11,1);phi2d0];
qd=EEE\RRHS;   
end

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);
 
% Data Storage Arrays
Q=zeros(ngc,10);
Qd=zeros(ngc,10);
Qdd=zeros(ngc,10);
Qddd=zeros(ngc,10);
Qh=zeros(ngc,10);
Qdh=zeros(ngc,10);
Qddh=zeros(ngc,10);
LLam=zeros(nc,10);
LLamd=zeros(nc,10);

%Position Correction
z=zeros(ngc,1); %z plays the role of deltaq
nu=zeros(nc,1);
Pz=[eye(ngc),zeros(ngc,nc)];
Pnu=[zeros(nc,ngc),eye(nc)];
err=intol+1;
i=0;
while err > intol
Phiq=PhiqEval(0,q+z,PJDT,par);
Resid=[z+Phiq'*nu;PhiEval(0,q+z,PJDT,par)];
I=eye(ngc);
JJ=[I+P4Eval(0,q+z,nu,PJDT,par),Phiq';Phiq,zeros(nc,nc)];
w=-JJ\Resid;
z=z+Pz*w;
nu=nu+Pnu*w;
err=norm(Resid);
i=i+1;
end
q=q+z;

%Velocity Correction
Phiq=PhiqEval(0,q,PJDT,par);
[Pst,Pstt,Pstq,Psttq]=P5Eval(0,q,par);
mu=(Phiq*Phiq')\(Phiq*qd+Pst);
delqd=-Phiq'*mu;
qd=qd+delqd;

Q(:,1)=q;
Qd(:,1)=qd;

%Calculate Initial Acceleration and Lagrange Multipliers
M=MEval(PMDT,par);
Phiq=PhiqEval(0,q,PJDT,par);
QA=QAEval(0,q,qd,PMDT,PTSDAT,PRSDAT,par);
Gam=GamEval(0,q,qd,PJDT,par);
EE=[M,Phiq';Phiq,zeros(nc,nc)];
CondEE=cond(EE);
RHS=[QA;-Gam];
x=EE\RHS;
Pqdd=[eye(ngc),zeros(ngc,nc)];
PLam=[zeros(nc,ngc),eye(nc)];
qdd=Pqdd*x;
Lam=PLam*x;

Qdd(:,1)=qdd;
LLam(:,1)=Lam;

%Initialize Data For Integration
n=1;
t(1)=0;          %Reparameterize if Cr>1
Jiter=0;
Jterurpt(1)=0;
JCondrpt(1)=0;
R1Normrpt(1)=R1nmax+1;
ECondrpt(1)=MaxECond+1;
Econd(1)=1;
VelConstrNorm(1)=0;
AccConstrNorm(1)=0;
corvel=0;
corpos=0;
coracc=0;
nch=1;

%Integration
while t(n)<=tfinal
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1); 

% Integration
if integ==1     %Implicit Index 1    
[q,qd,qdd,qdddnm,Lam,Lamdnm,R1n,Jiter,JCond,h,nch,Err]=...
    ImplicitIndex1(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
    PMDT,PTSDAT,PRSDAT,PJDT,par,alpha,nch);

R1Normrpt(n)=R1n;
Jiterrpt(n)=Jiter;
JCondrpt(n)=JCond;
Errrpt(n)=Err;
hrpt(n)=h;
end

if integ==2     %Implicit Index 2 
[q,qd,qdd,qdddnm,Lam,Lamdnm,R1n,Jiter,JCond,h,nch,Err]=...
    ImplicitIndex2(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
    PMDT,PTSDAT,PRSDAT,PJDT,par,alpha,nch);
    
R1Normrpt(n)=R1n;
Jiterrpt(n)=Jiter;
JCondrpt(n)=JCond;
Errrpt(n)=Err;
hrpt(n)=h;
end

if integ==3     %Implicit Index 3 
[q,qd,qdd,qdddnm,Lam,Lamdnm,R1n,Jiter,JCond,h,nch,Errn,hopt]=...
ImplicitIndex3(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
PMDT,PTSDAT,PRSDAT,PJDT,par,alpha,nch);
 
R1Normrpt(n)=R1n;
Jiterrpt(n)=Jiter;
JCondrpt(n)=JCond;
Errrpt(n)=Errn;
hrpt(n)=h;
hoptrpt(n)=hopt;
end

if integ==4     %Explicit RKFN45  
[q,qd,qdd,ECond,Err,h,nch]=ExplicitRKFN45(n,tn,Q,Qd,h,hmax,par,...
    PMDT,PTSDAT,PRSDAT,PJDT,nch);
qdddnm=zeros(ngc,1);
Lamdnm=zeros(nc,1);
ECondrpt(n)=ECond;
Errrpt(n)=Err;
hrpt(n)=h;
end

if integ==5     %Explicit Nystrom4 
[q,qd,qdd,Lam,ECond]=ExplicitNystrom4(n,tn,Q,Qd,...
    h,PMDT,PTSDAT,PRSDAT,PJDT,par);
qdddnm=zeros(ngc,1);
Lamdnm=zeros(nc,1);
ECondrpt(n)=ECond;    
end

if constrcor==1
%Corrections if position, velocity, or acceleration constraint 
%errors exceed tolerances

%Position Correction
if norm(PhiEval(tn,q,PJDT,par))>PosConstrMax
z=zeros(ngc,1); %z plays the role of deltaq
nu=zeros(nc,1);
Pz=[eye(ngc),zeros(ngc,nc)];
Pnu=[zeros(nc,ngc),eye(nc)];
err=intol+1;
i=0;
while err > intol
Phiq=PhiqEval(tn,q+z,PJDT,par);
Resid=[z+Phiq'*nu;PhiEval(tn,q+z,PJDT,par)];
I=eye(ngc);
JJ=[I+P4Eval(tn,q+z,nu,PJDT,par),Phiq';Phiq,zeros(nc,nc)];
w=-JJ\Resid;
z=z+Pz*w;
nu=nu+Pnu*w;
err=norm(Resid);
i=i+1;
end
q=q+z;
corpos=corpos+1;
corposrpt(n)=corpos;
corpositer(corpos)=i;
end

%Velocity Correction
Phiq=PhiqEval(tn,q,PJDT,par);
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
if norm(Phiq*qd+Pst)>VelConstrMax 
mu=(Phiq*Phiq')\(Phiq*qd+Pst);
delqd=-Phiq'*mu;
qd=qd+delqd;
corvel=corvel+1;
corvelrpt(n)=corvel;
end

%Acceleration Correction
Gam=GamEval(tn,q,qd,PJDT,par);
if norm(Phiq*qdd+Gam)>AccConstrMax
M=MEval(PMDT,par);
QA=QAEval(0,q0,qd0,PMDT,PTSDAT,PRSDAT,par);
EEE=[M,Phiq';Phiq,zeros(nc,nc)];
CondEE=cond(EEE);
RHS=[QA;-Gam];
x=EEE\RHS;
Pqdd=[eye(ngc),zeros(ngc,nc)];
PLam=[zeros(nc,ngc),eye(nc)];
qdd=Pqdd*x;
Lam=PLam*x;
coracc=coracc+1;
coraccrpt(n)=coracc;
end

end
%End Corrections


%Record Solution
Q(:,n)=q;
Qd(:,n)=qd;
Qdd(:,n)=qdd;
Qddd(:,n-1)=qdddnm;
qdddNorm(n)=norm(qdddnm);
LLam(:,n)=Lam;
LLamd(:,n-1)=Lamdnm;
LamdNorm(n)=norm(Lamdnm);

%Calculate Total Energy
%Kinetic Energy
M=MEval(PMDT,par);
KE(n)=0.5*qd'*M*qd;
%Potential Energy of Gravitational Force
PEg=0;
i=1;
while i<=nb
m=PMDT(1,i);
PEg=PEg+m*g*q(2+3*(i-1));
i=i+1;
end

%Potential Energy of Spring Force
PEs=0;
if NTSDA>0
T=1;
while T<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=PTSDATPart(PTSDAT,T);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
el=sqrt(d12'*d12);
%K=K*(1-sign(q(1)))/2;      %Unilateral spring constant, Double pendulum
PEs=PEs+0.5*K*(el-el0)^2;
T=T+1;
end

end

if NRSDA>0
R=1;
while R<=NRSDA
[i,j,K,C,phi0,T]=PRSDATPart(PRSDAT,R);
[r1,ph1]=qPart(q,i);
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
PEs=PEs+0.5*K*(ph2-ph1-phi0)^2;
R=R+1;
end
end

PE(n)=PEg+PEs;

TE(n)=KE(n)+PE(n);

%Data of Interest (Enter for each application)

if app==1
phi1(n)=q(3);
phi2(n)=q(6);
phi2d(n)=qd(6);
x4(n)=q(10);
x4d(n)=qd(10);
x4dd(n)=qdd(10);
Lam1rpt(n)=Lam(1);
end

if app==2
y1(n)=q(2);
y2(n)=q(5);
phi1(n)=q(3);
phi2(n)=q(6);    
end

if app==3
phi1(n)=q(3);
phi1d(n)=qd(3);
x2(n)=q(4); 
x2d(n)=qd(4);
x2dd(n)=qdd(4);
end

if app==4
phi1(n)=q(3);
x2(n)=q(4);    
end

if app==5
x1(n)=q(1); 
y1(n)=q(2); 
phi1(n)=q(3);   
end

if app==6
x1(n)=q(1); 
y1(n)=q(2); 
phi1(n)=q(3);   
end

if app==7
x1(n)=q(1); 
y1(n)=q(2); 
phi1(n)=q(3); 
x2(n)=q(4); 
y2(n)=q(5); 
phi2(n)=q(6);
phi2m1(n)=phi2(n)-phi1(n); 
end

if app==8
x1(n)=q(1); 
y1(n)=q(2); 
phi1(n)=q(3); 
x2(n)=q(4); 
y2(n)=q(5); 
phi2(n)=q(6);
phi2m1(n)=phi2(n)-phi1(n);
end

if app==9 
x1del(n)=q(1)-0.2; 
x2del(n)=q(4)-0.4;
x3del(n)=q(7)-0.6;
x4del(n)=q(10)-0.8;
x5del(n)=q(13)-1;

x1d(n)=qd(1); 
x2d(n)=qd(4);
x3d(n)=qd(7);
x4d(n)=qd(10);
x5d(n)=qd(13);

x1dd(n)=qd(1); 
x2dd(n)=qd(4);
x3dd(n)=qd(7);
x4dd(n)=qd(10);
x5dd(n)=qd(13);
end

if app==10 
x1del(n)=q(1)-0.1; 
x2del(n)=q(4)-0.2;
x3del(n)=q(7)-0.3;
x4del(n)=q(10)-0.4;
x5del(n)=q(13)-0.5;
x6del(n)=q(16)-0.6; 
x7del(n)=q(19)-0.7;
x8del(n)=q(22)-0.8;
x9del(n)=q(25)-0.9;
x10del(n)=q(28)-1;

x1d(n)=qd(1); 
x2d(n)=qd(4);
x3d(n)=qd(7);
x4d(n)=qd(10);
x5d(n)=qd(13);
x6d(n)=qd(16); 
x7d(n)=qd(19);
x8d(n)=qd(22);
x9d(n)=qd(25);
x10d(n)=qd(28);

x1dd(n)=qd(1); 
x2dd(n)=qd(4);
x3dd(n)=qd(7);
x4dd(n)=qd(10);
x5dd(n)=qd(13);
x6dd(n)=qd(16); 
x7dd(n)=qd(19);
x8dd(n)=qd(22);
x9dd(n)=qd(25);
x10dd(n)=qd(28);
end

%Calculate constraint error
Phi=PhiEval(tn,q,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
Gam=GamEval(tn,q,qd,PJDT,par);
PosConstrNorm(n)=norm(Phi);
VelConstrNorm(n)=norm(Phiq*qd);
AccConstrNorm(n)=norm(Phiq*qdd+Gam);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
EOMNorm(n)=norm(M*qdd+Phiq'*Lam-QA);

end



