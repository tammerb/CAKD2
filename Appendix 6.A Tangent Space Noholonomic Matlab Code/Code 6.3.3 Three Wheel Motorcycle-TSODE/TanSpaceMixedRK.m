
%Tangent Space Formulation-Tyicycle


utol=0.0001;     %Tolerance in solving for u
Btol=0.0001;    %Convergence criteria in B iteration
Htol=0.0001;    %Convergence criteria in H iteration
Maxv=20;         %Limit on magnitude of v
MaxRMcond=100;     %Limit on magnitude of reduced mass matrix condition number
Maxuiter=4;     %Limit on u iterations
MaxBiter=4;     %Limit on B iterations
MaxHiter=4;     %Limit on H iterations
h=0.005;       %Initial step Size
h2=h^2;
hmax=0.1;       %Maximum step size
tfinal=10;

integ=2;        %integ=1, RK4; integ=2, RKF45

nq=10;
nh=6;
nd=2;
nv=nq-nh;
nw=nq-nh-nd;
nu=nh;
nx=nh+nd;

m=200;
g=9.8;
k=15000;   %Spring Constant
c=1500;   %Damping Constant
F=0;   %Appllied force/torque
phi=0.2;        %phi=  0.2,Standard;  pi/2, chopper
amp=-0.004;    %Thet=amp*sin(om*t), 0<om*t<2*pi; zero after
om=1;
mode=1;     %mode=1, lane change; mode=2, step steer; mode=3, continuous sin

sf=(0.6+(m*g)/(2*k))/cos(phi);

%Enter all parameters to be used
par=[nq;nh;nd;nv;nu;nw;nx;m;g;F;k;c;phi;amp;om;sf;utol;Btol;Htol;mode;...
    integ;hmax];


Jv=[30;15;30];
p1=[cos(phi/2);sin(phi/2);0;0];
%Data to be passed by argument:
      %J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa
J=diag(Jv);
dpP0=[0;-0.6;0];
dpP1=[0.5;-0.6;0];
dpP2=[-0.5;-0.6;0];
bp=[0;0.6;0.6];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
P=[0,-1;1,0];
A1=ATran(p1);
apppsa=[0,0;eye(2)];
atpppsa=[0,0;P];

%q=[r;p;a;s]

% Data Storage Arrays
Q=zeros(nq,10);
Qd=zeros(nq,10);
Qdd=zeros(nq,10);
vrpt=zeros(nv,10);
vdrpt=zeros(nv,10);
vddrpt=zeros(nv,10);
urpt=zeros(nu,10);
udrpt=zeros(nu,10);
uddrpt=zeros(nu,10);
wrpt=zeros(nw,10);
wdrpt=zeros(nw,10);
wddrpt=zeros(nw,10);
xrpt=zeros(nx,10);
xdrpt=zeros(nx,10);
xddrpt=zeros(nx,10);
yf=zeros(nv+nw,10);
QArpt=zeros(nq,10);
atrpt=zeros(3,10);
Lam=zeros(nx,10);
QC3=zeros(nq,10)


% Initial Conditions
t0=0;
x0=0;
y0=0;
z0=0.304;
r0=[x0;y0;z0];
p0=[1;0;0;0];
a0=[-0.062;-0.294];
s0=0.621
q0=[r0;p0;a0;s0];


xd0=0;
yd0=15;
zd0=0;
sd0=0;
pd0=[0;0;0;0];
ad0=[0;0];
sd0=0;
qd0=[xd0;yd0;zd0;pd0;ad0;sd0];

Q(:,1)=q0;
Qd(:,1)=qd0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);
RMcond(1)=0;

%Start Integration Process
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
RMcond(1)=MaxRMcond+1;
uiterrpt(1)=Maxuiter+1;         
Biterrpt(1)=MaxBiter+1; 
Hiterrpt(1)=MaxHiter+1;
jRepar=0;       %Counter for Reparameterization
energylost(1)=0;

while t(n)<tfinal;

n=n+1
t(n)=t(n-1)+h;
tn=t(n);

%Criteria to Enter Parameterization 

if Biterrpt(n-1)>MaxBiter
    Cr=2;
end

if Hiterrpt(n-1)>MaxHiter
    Cr=2;
end

if uiterrpt(n-1)>Maxuiter
    Cr=2;
end

if vnorm(n-1)>Maxv
    Cr=2;
end

if RMcond(n-1)>MaxRMcond
    Cr=Cr+2;
end

if Cr>1     %Criteria for Parameterization
 
Crrpt(n)=Cr;
    
% Parameterization
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);
tnm=t(n-1);

[q0,qd0,V,U,W,X,B,H,jRepar]=Param(tnm,qnm,qdnm,jRepar,par,J,...
    dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);


%Enter v nvx1 zeros and u nux1 zeros
v=[0;0;0;0];   
u=[0;0;0;0;0;0];
w=[0;0];
x=[0;0;0;0;0;0;0;0];
vrpt(:,n-1)=v;
urpt(:,n-1)=u;
wrpt(:,n-1)=w;
xrpt(:,n-1)=x;
yf(:,n-1)=[v;w]; 


jReparrpt(n)=jRepar;
end


Cr=0

% Integration

t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1)
ue=urpt(:,n-1);
ynm=yf(:,n-1);



if integ==1
[yn,ydn,RM]=RungeKutta4FirstOrderIntegrate(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,h,h2,utol,Btol,Htol,par,...
    J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
end

if integ==2
[yn,ydn,RM,h,h2,err]=RKF45FirstOrderIntegrate(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,h,h2,utol,Btol,Htol,par,...
    J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
error(n)=err;
hrpt(n)=h;
end

RMcond(n)=cond(RM);

yf(:,n)=yn;
[v,w] = yPart(yn);
vn=v;
wn=w;
vrpt(:,n)=vn;
vnorm(n)=norm(vn);
wrpt(:,n)=wn;
[vd,wd] = ydPart(ydn);
vdn=vd;
wdn=wd;
vdrpt(:,n)=vdn;
vdnorm(n)=norm(vdn);
wdrpt(:,n)=wdn;


%Evaluate/update u, q, B, and H
[u,uiter]=usolv(tn,ue,vn,q0,V,U,B,utol,par,J,dpP0,dpP1,dpP2,...
    bp,ux,uy,uz,P,A1,apppsa,atpppsa);
un=u;
urpt(:,n)=un;
unorm(n)=norm(un);
uiterrpt(n)=uiter;
qn=q0+V*vn-U*un;
Q(:,n)=qn;
qnorm(n)=norm(qn);
[B,Biter]=Bcorr(tn,qn,B,U,Btol,par,J,dpP0,dpP1,dpP2,bp,...
    ux,uy,uz,P,A1,apppsa,atpppsa);
Biterrpt(n)=Biter;
C=Ceval(tn,qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
[H,Hiter]=Hcorr(H,X,C,Htol,par,J,dpP0,dpP1,dpP2,bp,...
    ux,uy,uz,P,A1,apppsa,atpppsa);
Hiterrpt(n)=Hiter;

%Evaluate qd and qdd

Nu=Nueval(tn,qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,...
    apppsa,atpppsa);
D=(eye(nq)-X*H*C)*W;
qdn=D*wn+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Qd(:,n)=qdn;
qdnorm(n)=norm(qdn);
Phiq=P1(tn,qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
Pt=-[Nu(1);Nu(2);Nu(3);Nu(4);Nu(5);Nu(6)];
udn=B*Phiq*V*vdn+B*Pt;
udrpt(:,n)=udn;
udnorm(n)=norm(udn);
Gam=Gameval(tn,qn,qdn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,...
    uz,P,A1,apppsa,atpppsa);
qddn=D*wdn-X*H*Gam;
Qdd(:,n)=qddn;
qddnorm(n)=norm(qddn);


%Data for further calculation
[r,p,a,s]=qPart(qn);
rn=r;
pn=p;
an=a;
sn=s;
[rd,pd,ad,sd]=qdPart(qdn);
rdn=rd;
pdn=pd;

rx(n)=rn(1);
ry(n)=rn(2);
rz(n)=rn(3);
srpt(n)=sn;

QA=QAEval(tn,qn,qdn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,...
    P,A1,apppsa,atpppsa)
QArpt(:,n)=QA
atrpt(:,n)=ATran(pn)*A1*[0;P*an]


%Kinetic and Total Energy
M=Meval(qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
KE(n)=0.5*qdn'*M*qdn;
TE(n)=KE(n)+0.5*k*(sn-sf)^2+m*g*rn(3);

%Constraint reaction forces
M=Meval(qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
S=Seval(qn,qdn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
QA=QAEval(tn,qn,qdn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,...
    P,A1,apppsa,atpppsa);
C=Ceval(tn,qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
Lambda=H'*X'*(M*qddn+S-QA);
Lam(:,n)=Lambda;
FNLR(n)=Lambda(2);
FNRR(n)=Lambda(1);
FNF(n)=Lambda(3);
FXpR(n)=Lambda(7);
FXpppF(n)=Lambda(8);
muf(n)=Lambda(8)/Lambda(3);
mur(n)=Lambda(7)/Lambda(2);


%Constraint error report
Phi=P0(tn,qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
phiErr(n)=norm(Phi);
C=Ceval(tn,qn,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
VelErr(n)=norm(C*qdn-Nu);
AccErr(n)=norm(C*qddn+Gam);



end




   