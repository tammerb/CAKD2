
%Tangent Space MixHolNonhol-Three Wheeled Transporter


Maxv=100;         %Limit on magnitude of v
Maxw=1;         %Limit on magnitude of w
MaxRMcond=10;     %Limit on magnitude of reduced mass matrix condition number

h=0.001;       %Step Size
hmax=0.01;
h2=h^2;
tfinal=10;

integ=1;        %integ=1, RK4; integ=2, RKF45

m=100;  %Moment of inertia also 10
g=9.8;
FL=0;
FR=0;

nq=3;
nh=0;
nd=1;
nv=nq-nh;
nu=nh;
nw=nq-nh-nd;
nx=nh+nd;


%Enter all parameters to be used
par=[nq;nh;nd;nv;nu;nw;nx;m;g;FL;FR;integ];

[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);

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
Lam=zeros(nh+nd,10);
QC=zeros(nq,10);
QC1=zeros(nq,10);
QC2=zeros(nq,10);

% Initial Conditions;
q0=[0;0;0];
qd0=[0;5;0];

Q(:,1)=q0;
Qd(:,1)=qd0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);


%Start Integration Process
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
wnorm(1)=Maxw+1;
RMcond(1)=MaxRMcond+1;
jRepar=0;       %Counter for Reparameterization

energylost(1)=0;

while t(n)<tfinal

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

%Criteria to Enter Parameterization 

if vnorm(n-1)>Maxv
    Cr=2;
end

if wnorm(n-1)>Maxw
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

[q0,qd0,V,W,X,jRepar]=Param(tnm,qnm,qdnm,jRepar,par);


%Enter v nvx1 zeros and u nux1 zeros
v=[0;0;0];   
w=[0;0];
x=0;
vrpt(:,n-1)=v;
wrpt(:,n-1)=w;
xrpt(:,n-1)=x;
yf(:,n-1)=[v;w]; 


jReparrpt(n)=jRepar;
end


Cr=0;

% Integration

tnm=t(n-1);
ue=0;
ynm=yf(:,n-1);

if integ==1
[yn,ydn,RM]=RungeKutta4FirstOrderIntegrate(tnm,ynm,...
    ue,V,W,X,q0,qd0,h,h2,par);
end

if integ==2
   
[yn,ydn,RM,h,h2,err]=RKF45FirstOrderIntegrate(tnm,ynm,...
    ue,V,W,X,q0,qd0,h,h2,hmax,par);

error(n)=err;
hrpt(n)=h;
end


RMcond(n)=cond(RM);

ynorm(n)=norm(yn);

yf(:,n)=yn;
[v,w] = yPart(yn);
vn=v;
wn=w;
vrpt(:,n)=vn;
vnorm(n)=norm(vn);
wnorm(n)=norm(wn);
wrpt(:,n)=wn;
[vd,wd] = ydPart(ydn);
vdn=vd;
wdn=wd;
vdrpt(:,n)=vdn;
vdnorm(n)=norm(vdn);
wdrpt(:,n)=wdn;


%Evaluate/update  q


qn=q0+V*vn;
Q(:,n)=qn;
qnorm(n)=norm(qn);

%Evaluate qd and qdd

C=Ceval(tn,qn,par);
a=C*X;
D=(eye(nq)-(1/a)*X*C)*W;
qdn=D*wn+(eye(nq)-(1/a)*X*C)*qd0;
Qd(:,n)=qdn;
qdnorm(n)=norm(qdn);

Gam=Gameval(tn,qn,qdn,par);
qddn=D*wdn+(Gam/a)*X;
Qdd(:,n)=qddn;
qddnorm(n)=norm(qddn);

%Kinetic and Total Energy
M=Meval(qn,par);
KE(n)=0.5*qdn'*M*qdn;

%Quantities of Interest
[r,phi]=qPart(qn);
rn=r;
phin=phi;
[rd,phid]=qdPart(qdn);
rdn=rd;
phidn=phid;

rx(n)=r(1);
ry(n)=r(2);


%Constraint error report

C=Ceval(tn,qn,par);
VelErr(n)=norm(C*qdn);
AccErr(n)=norm(C*qddn-Gam);



end




   