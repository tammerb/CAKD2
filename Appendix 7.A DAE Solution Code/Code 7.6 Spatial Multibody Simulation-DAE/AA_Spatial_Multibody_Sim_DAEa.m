%AA Spatial Multibody Sim-DAE

intol=10^-6;     %Tolerance in solving discretized equations of motion
Atol=10^-5;     %Absolute error tolerance for variable step metnods
MaxIntiter=8;      %Limit on number of integration iterations
MaxJcond=200;   %Limit on magnitude of Jcond
R1nmax=15;     %Limit on residual in integrationt
MaxECond=10;  %Limit on magnitude of ECond
PosConstrMax=10^-4;  %Limit on position constraint error
VelConstrMax=10^-4;  %Limit on velocity constraint error
AccConstrMax=10^20;  %Limit on acceleration constraint error

h0=0.001;       %Initial time step
h=h0;
hmax=0.001;
hvar=1;     %hvar=1, variable h;hvar=2, constant h

constrcor=1;         %constrcor=1, correct; constrcor=2, no correct

tfinal=100;

integ=1;    %integ=1-ImplicitIndex1; alpha=0-Trapezoidal; alpha<0-HHt
            %integ=2-ImplicitIndex2; alpha=0-Trapezoidal; alpha<0-HHt
            %integ=3-ImplicitIndex3; alpha=0-Trapezoidal; alpha<0-HHt
            %integ=4-ExplicitNystrom4
            %integ=5-ExplicitRKFN45
            
alpha=-1/3;    %Enter -1/3<=alpha<0 for HHT            
                %alpha not used in ExplicitNystrom4 and ExplicitRKFN45          
                       
g=9.8;

%Application Data
    %app=1, Pendulum, Spherical to Ground
    %app=2, Top, Spherical to Ground
    %app=3, One Body Pendulum, Dist. to Ground
    %app=4  Double Pendulum
    %app=5, One Body Cylindrical with Spring
    %app=6, Spatial Slider-Crank
    
app=3;

[nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q,qd]=AppData(app);

par=[nb;ngc;nh;nc;g;intol;Atol;h0;hvar;NTSDA];

%Initial condition calculation, if required

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

% Data Storage Arrays
Q=zeros(ngc,10);
Qd=zeros(ngc,10);
Qdd=zeros(ngc,10);
LLam=zeros(nc,10);
ConPhi=zeros(nc,10);

%Initial position correction
z=zeros(ngc,1);
nu=zeros(nc,1);
Pz=[eye(ngc),zeros(ngc,nc)];
Pnu=[zeros(nc,ngc),eye(nc)];
err=intol+1;
i=0;
while err > intol
Phiq=PhiqEval(0,q+z,SJDT,par);
Resid=[z+Phiq'*nu;PhiEval(0,q+z,SJDT,par)];
I=eye(ngc);
JJ=[I+P4Eval(0,q+z,nu,SJDT,par),Phiq';Phiq,zeros(nc,nc)];
w=-JJ\Resid;
z=z+Pz*w;
nu=nu+Pnu*w;
err=norm(Resid);
i=i+1;
end
q=q+z;

%Initial velocity correction
Phiq=PhiqEval(0,q,SJDT,par);
mu=(Phiq*Phiq')\Phiq*qd;
delqd=-Phiq'*mu;
qd=qd+delqd;

Q(:,1)=q;
Qd(:,1)=qd;

%Calculate Initial Acceleration and Lagrange Multipliers
M=MEval(q,SMDT,par);
Phiq=PhiqEval(0,q,SJDT,par);
QA=QAEval(0,q,qd,SMDT,STSDAT,par);
S=SEval(q,qd,SMDT,par);
Gam=GamEval(0,q,qd,SJDT,par);
EE=[M,Phiq';Phiq,zeros(nc,nc)];
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
Jiterrpt(1)=0;
Jiter=0;
JCondrpt(1)=0;
R1Normrpt(1)=R1nmax+1;
ECondrpt(1)=MaxECond+1;
Econd(1)=1;
PosConstrNorm(1)=0;
VelConstrNorm(1)=0;
AccConstrNorm(1)=0;
corvel=0;
corpos=0;
nch=1;

%Integration

while t(n)<tfinal

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

% Integration

if integ==1
[q,qd,qdd,Lam,R1n,Jiter,JCond,h,nch,Err]=...
    ImplicitIndex1(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
    SMDT,STSDAT,SJDT,par,alpha,nch);

R1Normrpt(n)=R1n;
Jiterrpt(n)=Jiter;
JCondrpt(n)=JCond;
Errrpt(n)=Err;
hrpt(n)=h;

end

if integ==2
[q,qd,qdd,Lam,R1n,Jiter,JCond,h,nch,Err]=...
    ImplicitIndex2(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
    SMDT,STSDAT,SJDT,par,alpha,nch);

R1Normrpt(n)=R1n;
Jiterrpt(n)=Jiter;
JCondrpt(n)=JCond;
Errrpt(n)=Err;
hrpt(n)=h;
end

if integ==3
[q,qd,qdd,Lam,R1n,Jiter,JCond,h,nch,Err,hopt]=...
    ImplicitIndex3(n,tn,Q,Qd,Qdd,LLam,h,hmax,...
    SMDT,STSDAT,SJDT,par,alpha,nch);

R1Normrpt(n)=R1n;
Jiterrpt(n)=Jiter;
JCondrpt(n)=JCond;
Errrpt(n)=Err;
hrpt(n)=h;
end

if integ==4   
[q,qd,qdd,Lam,ECond]=ExplicitNystrom4(n,tn,Q,Qd,h,...
    SMDT,STSDAT,SJDT,par);
ECondrpt(n)=ECond;    
end

if integ==5
[q,qd,qdd,ECond,h,nch]=ExplicitRKFN45(n,tn,Q,Qd,h,hmax,par,...
    SMDT,STSDAT,SJDT,nch);
ECondrpt(n)=ECond;
hrpt(n)=h;
end

%Corrections if velocity or position errors exceed tolerances
if constrcor==1

Phiq=PhiqEval(tn,q,SJDT,par);
if norm(Phiq*qd)>VelConstrMax 
mu=(Phiq*Phiq')\Phiq*qd;
delqd=-Phiq'*mu;
qd=qd+delqd;
corvel=corvel+1;
corvelrpt(n)=corvel;
end

if norm(PhiEval(tn,q,SJDT,par))>PosConstrMax
z=zeros(ngc,1);
nu=zeros(nc,1);
Pz=[eye(ngc),zeros(ngc,nc)];
Pnu=[zeros(nc,ngc),eye(nc)];
err=intol+1;
i=0;
while err > intol
Phiq=PhiqEval(tn,q+z,SJDT,par);
Resid=[z+Phiq'*nu;PhiEval(tn,q+z,SJDT,par)];
I=eye(ngc);
JJ=[I+P4Eval(tn,q+z,nu,SJDT,par),Phiq';Phiq,zeros(nc,nc)];
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

end
%End Corrections

%Record Solution
Q(:,n)=q;
Qd(:,n)=qd;
Qdd(:,n)=qdd;
LLam(:,n)=Lam;

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

if app==1
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
end

if app==2
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
end

if app==3
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
omeg1norm(n)=norm(2*EEval(p1)*p1d);
end

if app==4
x2(n)=q(8);
y2(n)=q(9);
z2(n)=q(10);
x2d(n)=qd(8);
y2d(n)=qd(9);
z2d(n)=qd(10);
x2dd(n)=qdd(8);
y2dd(n)=qdd(9);
z2dd(n)=qdd(10);
end

if app==5
x1(n)=q(1);
y1(n)=q(2);
z1(n)=q(3);
end

if app==6
z2(n)=q(10);
z2d(n)=qd(10);
z2dd(n)=qdd(10);
[r1,p1]=qPart(q,1);
[r1d,p1d]=qPart(qd,1);
[r1dd,p1dd]=qPart(qdd,1);
omeg1=(2*EEval(p1)*p1d);
omeg1z(n)=omeg1'*uz;
omeg1zd(n)=(2*EEval(p1)*p1dd)'*uz;
end

%Calculate constraint error
Phi=PhiEval(tn,q,SJDT,par);
Phiq=PhiqEval(tn,q,SJDT,par);
Gam=GamEval(tn,q,qd,SJDT,par);
PosConstrNorm(n)=norm(Phi);
VelConstrNorm(n)=norm(Phiq*qd);
AccConstrNorm(n)=norm(Phiq*qdd+Gam);

end


