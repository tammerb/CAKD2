
%Tangent Space Formulation-Rolling Ellipsoid on moving surface


utol=0.0001;     %Tolerance in solving for u
Btol=0.0001;    %Convergence criteria in B iteration
Htol=0.0001;    %Convergence criteria in H iteration
Maxv=1;         %Limit on magnitude of v
MaxRMcond=10;     %Limit on magnitude of reduced mass matrix condition number
Maxuiter=4;     %Limit on u iterations
MaxBiter=4;     %Limit on B iterations
MaxHiter=4;     %Limit on H iterations
h=0.001;        %Step Size
hmax=0.01;
h2=h^2;
tfinal=10;

Integ=1;        %Integ=1, RKF; Integ=2, RK4


nq=12;
nh=7;
nd=2;
nv=nq-nh;
nw=nq-nh-nd;
nu=nh;
nx=nh+nd;

alph=0.4;
beta=0.4;
delt=0.4;


Lv=[1/(alph^2);1/(beta^2);1/(delt^2)];
L=diag(Lv);
g=9.8;
%density=140 kg/m3
m=4*140*pi*alph*beta*delt/3;
Jprv=(m/5)*[beta^2+delt^2;alph^2+delt^2;alph^2+beta^2];
Jpr=diag(Jprv);
K=150;   %Spring Constant
F=0;   %Appllied force is f=[F*sin(t),F*cos(2*t);0]
amp=0;
om=3;
eps=0;

%Enter all parameters to be used
par=[nq;nh;nd;nv;nu;nw;nx;m;g;F;K;eps;amp;om];

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

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
Ftan=zeros(2,10);

% Initial Conditions
t0=0;
x0=0;
y0=0;
p0=[1;0;0;0];
q0=[x0;y0;delt;1;0;0;0;0;0;-delt;x0;y0];


xd0=0.3;
yd0=0.3;
Omegz0=0;

%uz'*2*Ebareval(p0)=Omegz0
ICcoef=[Ceval(t0,q0,par,L);1,zeros(1,11);0,1,zeros(1,10);...
    0,0,0,uz'*2*Ebareval(p0),zeros(1,5)];
ICcond=cond(ICcoef);
Nu=Nueval(t0,q0,par);
ICrhs=[Nu;xd0;yd0;Omegz0];
qd0=ICcoef\ICrhs;

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

while t(n)<tfinal;

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

%Criteria to Enter Parameterization 

if Biterrpt(n-1)>MaxBiter
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

[q0,qd0,V,U,W,X,B,H,jRepar]=Param(tnm,qnm,qdnm,jRepar,par,L);

%Enter v nvx1 zeros and u nux1 zeros
v=[0;0;0;0;0];   
u=[0;0;0;0;0;0;0];
w=[0;0;0];
x=[0;0;0;0;0;0;0;0;0];
vrpt(:,n-1)=v;
urpt(:,n-1)=u;
wrpt(:,n-1)=w;
xrpt(:,n-1)=x;
yf(:,n-1)=[v;w]; 


jReparrpt(n)=jRepar;
end


Cr=0;

% Integration

t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);
ue=urpt(:,n-1);
ynm=yf(:,n-1);

if Integ==1
    
[yn,ydn,RM,h,h2,err,k]=RKF45FirstOrderIntegrate(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,h,h2,utol,Btol,Htol,par,L,Jpr,hmax);

hrpt(n)=h;
error(n)=err;

end

if Integ==2
    
[yn,ydn,RM]=RungeKutta4FirstOrderIntegrate(tnm,ynm,...
    ue,V,U,W,X,B,H,q0,qd0,h,h2,utol,Btol,Htol,par,L,Jpr);

end

RMcond(n)=cond(RM);
%error(n)=err;
hrpt(n)=h;

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
[u,uiter]=usolv(tn,ue,vn,q0,V,U,B,utol,par,L);
un=u;
urpt(:,n)=un;
unorm(n)=norm(un);
uiterrpt(n)=uiter;
qn=q0+V*vn-U*un;
Q(:,n)=qn;
qnorm(n)=norm(qn);
[B,Biter]=Bcorr(tn,qn,B,U,Btol,par,L);
Biterrpt(n)=Biter;
C=Ceval(tn,qn,par,L);
[H,Hiter]= Hcorr(H,X,C,Htol,par);
Hiterrpt(n)=Hiter;

%Evaluate qd and qdd

Nu=Nueval(tn,qn,par);
D2=(eye(nq)-X*H*C)*W;
qdn=D2*wn+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Qd(:,n)=qdn;
qdnorm(n)=norm(qdn);
Phiq=Phiqeval(tn,qn,par,L);
gtd=amp*om*sin(om*tn);
Pt=[-gtd*uz;0;0;0;0];
udn=B*Phiq*V*vdn+B*Pt;
udrpt(:,n)=udn;
udnorm(n)=norm(udn);
Gam=Gameval(tn,qn,qdn,par,L);
qddn=D2*wdn-X*H*Gam;
Qdd(:,n)=qddn;
qddnorm(n)=norm(qddn);


%Data for further calculation
[r,p,apr,x1,y1]=qPart(qn);
rn=r;
pn=p;
aprn=apr;
x1n=x1;
y1n=y1;
[rd,pd,aprd,x1d,y1d]=qdPart(qdn);
rdn=rd;
pdn=pd;

%Angular Velocity and Momentum

En=Ebareval(pn);
Gn=Gbareval(pn);
Omega(:,n)=2*En*pdn;
Omegapr(:,n)=2*Gn*pdn;
Angmom=ATran(pn)*Jpr*2*Gn*pdn;
Angmom(:,n)=Angmom;
zAngmom(n)=uz'*Angmom(:,n);
Mom(:,n)=m*rdn;
Momnorm(n)=norm(m*rdn);
if n==2
    Angmom(:,1)=Angmom(:,2);
    zAngmom(1)=zAngmom(2);
    Mom(:,1)=Mom(:,2);
    Momnorm(1)=Momnorm(2);
end

%Kinetic and Total Energy
M=Meval(qn,par,Jpr);
KE(n)=0.5*qdn'*M*qdn;
TE(n)=KE(n)+0.5*K*r'*r+m*g*rn(3);

%Quantities of Interest
rx(n)=rn(1);
ry(n)=rn(2);
rz(n)=rn(3);
rPn=rn+ATran(pn)*aprn;
rP(:,n)=rPn;
rPx(n)=rPn(1);
rPy(n)=rPn(2);
    

%Constraint Forces
tn=t(n);
f=[F*sin(tn);F*cos(2*tn);-m*g];
FA=[f-K*(rn-uz)];
rddn=[qddn(1);qddn(2);qddn(3)];
FC=m*rddn-FA;
nor=[-2*eps*x1;-4*eps*y1;1];
g1=[1;0;2*eps*x1];
g2=[0;1;4*eps*y1];
FC1=FC'*g1/norm(g1);
FC2=FC'*g2/norm(g2);
Ftan(:,n)=[FC1;FC2];
Ft=sqrt(FC1^2+FC2^2);
Ftannorm(n)=Ft;
Fnormal(n)=FC'*nor/norm(nor);
mu(n)=Ft/Fnormal(n);

if n==2
    rx(1)=rx(2);
    ry(1)=ry(2);
    Ftan(1)=Ftan(2);
    Fnormal(1)=Fnormal(2);
    mu(1)=mu(2);
    
end

%Constraint error report

Phi=Phieval(tn,qn,par,L);
phiErr(n)=norm(Phi);
Phiq=Phiqeval(tn,qn,par,L);
C=Ceval(tn,qn,par,L);
VelErr(n)=norm(C*qdn-Nu);
AccErr(n)=norm(C*qddn+Gam);



end




   