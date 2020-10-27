
%Tangent Space MixHolNonhol-Disk Rolling

utol=10^-8;     %Tolerance in solving for u
Htol=10^-8;    %Convergence criteria in H iteration
Btol=10^-8;    %Convergence criteria in B iteration
intol=10^-8;    %Integration error tolerance
Maxv=0.7;        %Limit on magnitude of v
MaxSolIter=6;  %Limit on number of integration iterations
Maxuiter=4;     %Limit on u iterations
MaxBiter=4;     %Limit on B iterations
MaxHiter=4;     %Limit on H iterations
MaxRMcond=100;     %Limit on RM condition number

h=0.0001;       %Step Size
hmax=0.0001;       %Maximum Allowed Step Size
h2=h^2;
tfinal=10;

integ=6;         %integ=1,Trap; integ=2,RK4; integ=3,SDIRK54; integ=4,BE1
                 %integ=5,BE4; integ=6,BE2

errcontr=2;      %errcontr=1, use in SDIRKL54; errcontr=2, not use in SDIRKL54

nq=9;
nh=4;
nd=2;
nv=nq-nh;
nw=nq-nh-nd;
nu=nh;
nx=nh+nd;

g=9.6;
m=10;
Jprv=[5;2.5001;2.5001];
Jpr=diag(Jprv);

%Enter all parameters to be used
par=[nq;nh;nd;nv;nu;nw;nx;m;g;Jprv;utol;Btol;Htol;intol;h];

apsa=[0,0;eye(2)];
bpsa=[0,0;0,-1;1,0];
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

% Data Storage Arrays
Q=zeros(nq,10);
Qd=zeros(nq,10);
Qdd=zeros(nq,10);
Vv=zeros(nv,10);
Vvd=zeros(nv,10);
Vvdd=zeros(nv,10);
Uu=zeros(nu,10);
Uud=zeros(nu,10);
Uudd=zeros(nu,10);
Ww=zeros(nw,10);
Wwd=zeros(nw,10);
Wwdd=zeros(nw,10);
Xx=zeros(nx,10);
Xxd=zeros(nx,10);
Xxdd=zeros(nx,10);
LLam=zeros(nh+nd,10);
Kk1=zeros(nv+nw,10);
Kk2=zeros(nv+nw,10);
Kk3=zeros(nv+nw,10);
Kk4=zeros(nv+nw,10);
Kk5=zeros(nv+nw,10);
Omeg=zeros(3,10);
uxrpt=zeros(3,10);

% Initial Conditions
tn(1)=0;
r0=uz;
p0=[1;0;0;0];
q0=[r0;p0;0;-1];

yd0=0;
omegx0=-1;
Omeg0=[omegx0;0.01;0];

%Calculate Initial Velocities
ICcoef=[Ceval(tn(1),q0,par);...
    zeros(3,3),2*Gbareval(p0),zeros(3,2)];
[Nu,Nusq]=NuNusqeval(tn(1),q0,par);
ICrhs=[Nu;Omeg0];
qd0=ICcoef\ICrhs;

%Calculate Total Energy
M=Meval(q0,par);
KE(1)=0.5*qd0'*M*qd0;
TE(1)=KE(1)+m*g;

%Calculate Initial qdd and Lam
M=Meval(q0,par);
QA=QAEval(tn(1),q0,qd0,par);
S=Seval(q0,qd0,par);
C=Ceval(tn(1),q0,par);
Gam=Gameval(tn(1),q0,qd0,par);

AccCoef=[M,C';C,zeros(6,6)];
AccRHS=[QA+S;-Gam];
z=AccCoef\AccRHS;

qdd0=[z(1);z(2);z(3);z(4);z(5);z(6);z(7);z(8);z(9)];
Lam0=[z(10);z(11);z(12);z(13);z(14);z(15)];

Q(:,1)=q0;
Qd(:,1)=qd0;
Qdd(:,1)=qdd0;
LLam(:,1)=Lam0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);

%Start Integration Process
n=1;
t1(1)=0;
vnorm(1)=Maxv+1;
SolIter(1)=MaxSolIter+1;
uiterrpt(1)=Maxuiter+1;         
Biterrpt(1)=MaxBiter+1; 
Hiterrpt(1)=MaxHiter+1;
RMcond(1)=MaxRMcond+1;
jRepar=0;       %Counter for Reparameterization

while t1(n)<tfinal

n=n+1;
t1(n)=t1(n-1)+h;
tn=t1(n);

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

if integ==1
if SolIter(n-1)>MaxSolIter
    Cr=Cr+2;
end
end

if RMcond(n-1)>MaxRMcond
    Cr=2;
end

if Cr>1     %Criteria for Parameterization
 
Crrpt(n)=Cr;
    
% Parameterization
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);
qddnm=Qdd(:,n-1);
tnm=t1(n-1);

[q0,qd0,V,U,W,X,B,H,jRepar]=Param(tnm,qnm,qdnm,jRepar,par);

jReparrpt(n)=jRepar;
nRepar=n-1;

%Enter v, u, w, and x zeros
vnm=zeros(nv,1);   
unm=zeros(nu,1);
wnm=zeros(nw,1);
xnm=zeros(nx,1);
Vv(:,n-1)=vnm;
Uu(:,n-1)=unm;
Ww(:,n-1)=wnm;
Xx(:,n-1)=xnm;

%Evaluate vdnm and wdnm

vdnm=V'*qdnm;
wdnm=W'*qddnm;
Vvd(:,n-1)=vdnm;
Wwd(:,n-1)=wdnm;

%Jacobian2 Evaluation
Lamnm=LLam(:,n-1);

J2=J2Eval(tn,qnm,qdnm,qd0,wnm,wdnm,Lamnm,U,V,W,X,B,H,par);

end

Cr=0;


% Integration
tnm=t1(n-1);
tn=tnm+h;
ue=Uu(:,n-1);

if integ==1
    
[vn,wn,vdn,wdn,Lam,RM,jiter,R1n]=TrapInteg(n,tn,Vv,Vvd,Ww,Wwd,ue,LLam,...
    V,U,W,X,B,H,q0,qd0,J2,par,nRepar);

SolIter(n)=jiter;
R1nrpt(n)=R1n;
end

if integ==2

vnm=Vv(:,n-1);
wnm=Ww(:,n-1);
       
[vn,wn,vdn,wdn,Lam,RM]=RungeKutta4FirstOrderIntegrate(tnm,...
   vnm,wnm,ue,V,U,W,X,B,H,q0,qd0,par);

end

if integ==3
    
[vn,vdn,wn,wdn,Lam,jiter,RM,k1,k2,k3,k4,k5,RR,h,hopt]=...
    SDIRK54Integ(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,...
    Kk5,q0,qd0,V,U,W,X,B,H,h,hmax,utol,intol,Btol,Htol,par,J2,errcontr);

Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
hrpt(n)=h;
hoptrpt(n)=hopt;
end

if integ==4
    
[vn,vdn,wn,wdn,Lam,jiter,RM,k1,RR]=...
    BE1Integ(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,q0,qd0,...
    V,U,W,X,B,H,h,utol,intol,Btol,Htol,par,J2,errcontr);

Kk1(:,n)=k1;

end

if integ==5
    
[vn,vdn,wn,wdn,Lam,jiter,RM,k1,k2,k3,k4,RR]=...
    BE4Integ(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,...
    q0,qd0,V,U,W,X,B,H,h,hmax,utol,intol,Btol,Htol,par,J2,errcontr);

Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
end

if integ==6
    
[vn,vdn,wn,wdn,Lam,jiter,RM,k1,k2,RR]=...
    BE2Integ(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,Kk2,...
    q0,qd0,V,U,W,X,B,H,h,hmax,utol,intol,Btol,Htol,par,J2,errcontr);

Kk1(:,n)=k1;
Kk2(:,n)=k2;

end

RMcond(n)=cond(RM);

Vv(:,n)=vn;
vnorm(n)=norm(vn);
Ww(:,n)=wn;
wnorm(n)=norm(wn);
Vvd(:,n)=vdn;
vdnorm(n)=norm(vdn);
Wwd(:,n)=wdn;
wdnorm(n)=norm(wdn);
LLam(:,n)=Lam;
Lamnorm(n)=norm(Lam);


%Evaluate/update u, q, B, and H
[u,uiter]=usolv(tn,ue,vn,q0,V,U,B,par);
un=u;
Uu(:,n)=un;
unorm(n)=norm(un);
uiterrpt(n)=uiter;
qn=q0+V*vn-U*un;
Q(:,n)=qn;
qnorm(n)=norm(qn);
[B,Biter]=Bcorr(tn,qn,B,U,par);
Biterrpt(n)=Biter;
C=Ceval(tn,qn,par);
[H,Hiter]= Hcorr(H,X,C,par);
Hiterrpt(n)=Hiter;

%Evaluate qd and qdd

D=(eye(nq)-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(tn,qn,par);
qdn=D*wn+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Qd(:,n)=qdn;
qdnorm(n)=norm(qdn);
Phiq=Phiqeval(tn,qn,par);
udn=B*Phiq*V*vdn;
Uud(:,n)=udn;
udnorm(n)=norm(udn);
Gam=Gameval(tn,qn,qdn,par);
qddn=D*wdn-X*H*Gam;
Qdd(:,n)=qddn;
qddnorm(n)=norm(qddn);

%Data of Interest
rx(n)=qn(1);
ry(n)=qn(2);
rz(n)=qn(3);
rdx(n)=qdn(1);
rdy(n)=qdn(2);
rdz(n)=qdn(3);
[r,p,a]=qPart(qn);
[rd,pd,ad]=qdPart(qdn);
omeg=2*Ebareval(p)*pd;
Omeg(:,n)=omeg;
uxrpt(:,n)=ATran(p)'*ux;

%Constraint Forces
[r,p,a]=qPart(qn);
[rd,pd,ad]=qdPart(qdn);
FN(n)=abs(Lam(1));
FS(n)=norm([Lam(5);Lam(6)]);
EffFrictCoef(n)=FS(n)/FN(n);
VelP(n)=norm(rd+Bbareval(p,apsa*a)*pd);

%Total Energy
M=Meval(qn,par);
KE(n)=0.5*qdn'*M*qdn;
TE(n)=KE(n)+m*g*rz(n);


%Constraint error report
Phi=Phieval(tn,qn,par);
phiErr(n)=norm(Phi);
C=Ceval(tn,qn,par);
VelErr(n)=norm(C*qdn);
AccErr(n)=norm(C*qddn+Gam);

end

maxTE=max(TE);
minTE=min(TE);



   