
%Tangent Space MixHolNonhol-PlanarArtVen


utol=0.0000000001;     %Tolerance in solving for u
Btol=0.0000000001;    %Convergence criteria in B iteration
Htol=0.0000000001;    %Convergence criteria in H iteration
Maxv=10;        %Limit on magnitude of v
MaxSolIter=10;  %Limit on number of integration iterations
Maxuiter=4;     %Limit on u iterations
MaxBiter=4;     %Limit on B iterations
MaxHiter=4;     %Limit on H iterations
intol=0.0001;    %Integration error tolerance
h=0.001;     %Step Size
hmax=0.001;    %Maximum allowed step size
h2=h^2;
tfinal=10;

integ=3;         %integ=1, Trap; integ=2, RK4; integ=3, SDIRK54 

errcontr=2;      %errcontr=1, use in SDIRKL54; errcontr=2, not use in SDIRKL54


m=1500;  %Moment of inertia also 
g=9.6;
amp=0.017;  %thet=amp*sin(t)
om=1;
mode=1;      %mode=1, lane change; mode=2, step steer; mode=3, sine steer

nq=6;
nh=2;
nd=3;
nv=nq-nh;
nu=nh;
nw=nq-nh-nd;
nx=nh+nd;

%Enter all parameters to be used
par=[nq;nh;nd;nv;nu;nw;nx;m;g;amp;om;mode;integ;h;utol;Btol;Htol];

% Data Storage Arrays
Q=zeros(nq,10);
Qd=zeros(nq,10);
Qdd=zeros(nq,10);
Vv=zeros(nv,10);
Vvd=zeros(nv,10);
Uu=zeros(nu,10);
Uud=zeros(nu,10);
Ww=zeros(nw,10);
Wwd=zeros(nw,10);
Xx=zeros(nx,10);
Xxd=zeros(nx,10);
LLam=zeros(nh+nd,10);
Ff=zeros(nv+nw,10);
Kk1=zeros(nv+nw,10);Kk1=zeros(nv+nw,10);
Kk2=zeros(nv+nw,10);
Kk3=zeros(nv+nw,10);
Kk4=zeros(nv+nw,10);
Kk5=zeros(nv+nw,10);


% Initial Conditions;
q0=[6;0;0;0;0;0];
qd0=[15;0;0;15;0;0];

Q(:,1)=q0;
Qd(:,1)=qd0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);

%Kinetic Energy
M=Meval(q0,par);
KE(1)=0.5*qd0'*M*qd0;

%Initial qdd and Lam
t1(1)=0;
M=Meval(q0,par);
QA=QAEval(t1(1),q0,qd0,par);
S=Seval(q0,qd0,par);
C=Ceval(t1(1),q0,par);
Gam=Gameval(t1(1),q0,qd0,par);

AccCoef=[M,C';C,zeros(5,5)];
AccRHS=[QA+S;-Gam];
z=AccCoef\AccRHS;

qdd0=[z(1);z(2);z(3);z(4);z(5);z(6)];
Lam0=[z(7);z(8);z(9);z(10);z(11)];

Qdd(:,1)=qdd0;
LLam(:,1)=Lam0;

%Start Integration Process
n=1;
vnorm(1)=Maxv+1;
SolIter(1)=MaxSolIter+1;
uiterrpt(1)=Maxuiter+1;         
Biterrpt(1)=MaxBiter+1; 
Hiterrpt(1)=MaxHiter+1;
jRepar=0;       %Counter for Reparameterization

while t1(n)<tfinal;

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
vnm=[0;0;0;0];   
unm=[0;0];
wnm=0;
xnm=[0;0;0;0;0];
Vv(:,n-1)=vnm;
Uu(:,n-1)=unm;
Ww(:,n-1)=wnm;
Xx(:,n-1)=xnm;

%Jacobian2 Evaluation
vdnm=V'*qdnm;
wdnm=W'*qddnm;
Vvd(:,n-1)=vdnm;
Wwd(:,n-1)=wdnm;
Lam=LLam(:,n-1);

J2=J2Eval(tn,q0,qd0,wnm,wdnm,Lam,U,V,W,X,B,H,q0,qd0,par);

end

Cr=0;

% Integration

tnm=t1(n-1);
tn=tnm+h;
ue=Uu(:,n-1);

if integ==1
[vn,wn,vdn,wdn,Lam,RM,jiter,R1n]=TrapInteg(n,tn,Vv,Vvd,Ww,Wwd,ue,LLam,...
    V,U,W,X,B,H,q0,qd0,J2,h,h2,utol,Btol,Htol,intol,par,nRepar);
SolIter(n)=jiter;
R1nrpt(n)=R1n;

end

if integ==2
vnm=Vv(:,n-1);
wnm=Ww(:,n-1);
[vn,wn,vdn,wdn,Lam,RM,f]=RungeKutta4FirstOrderIntegrate(tnm,...
    vnm,wnm,ue,V,U,W,X,B,H,q0,qd0,par);
Ff(:,n)=f;
end
 if integ==3
     
 [vn,vdn,wn,wdn,Lam,jiter,RM,k1,k2,k3,k4,k5,RR,h]=...
    IntegTanSpSDIRK54PartNR(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,...
    Kk2,Kk3,Kk4,Kk5,q0,qd0,V,U,W,X,B,H,h,hmax,utol,intol,Btol,Htol,...
    par,J2,errcontr);

Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
hrpt(n)=h;
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


%Evaluate/update u, q, B, and H
[unm,uiter]=usolv(tn,ue,vn,q0,V,U,B,utol,par);
un=unm;
Uu(:,n)=un;
unorm(n)=norm(un);
uiterrpt(n)=uiter;
qn=q0+V*vn-U*un;
Q(:,n)=qn;
qnorm(n)=norm(qn);
[B,Biter]=Bcorr(tn,qn,B,U,Btol,par);
Biterrpt(n)=Biter;
C=Ceval(tn,qn,par);
[H,Hiter]= Hcorr(H,X,C,Htol,par);
Hiterrpt(n)=Hiter;

%Evaluate qd and qdd

D=(eye(nq)-X*H*C)*W;
qdn=D*wn+(eye(nq)-X*H*C)*qd0;
Qd(:,n)=qdn;
qdnorm(n)=norm(qdn);
Phiq=Phiqeval(tn,qn,par);
udn=B*Phiq*V*vdn;
udrpt(:,n)=udn;
udnorm(n)=norm(udn);
Gam=Gameval(tn,qn,qdn,par);
qddn=D*wdn-X*H*Gam;
Qdd(:,n)=qddn;
qddnorm(n)=norm(qddn);

%Kinetic Energy
M=Meval(qn,par);
KE(n)=0.5*qdn'*M*qdn;

%Quantities of Interest
[r1,phi1,r2,phi2]=qPart(qn);
[r1d,phi1d,r2d,phi2d]=qdPart(qdn);

r1x(n)=r1(1);
r1y(n)=r1(2);
r2x(n)=r2(1);
r2y(n)=r2(2);
delphi(n)=phi1-phi2;

%Constraint Forces
M=Meval(qn,par);
QA=QAEval(tn,qn,qdn,par);
Lambda=H'*X'*(M*qddn-QA);
Lam(:,n)=Lambda;
C=Ceval(tn,qn,par);
Qrev=C'*[Lambda(1);Lambda(2);0;0;0];
Qfw=C'*[0;0;Lambda(3);0;0];
Qiw=C'*[0;0;0;Lambda(4);0];
Qrw=C'*[0;0;0;0;Lambda(5)];
Frev(n)=norm([Qrev(1);Qrev(2)]);
Ffw(n)=norm([Qfw(1);Qfw(2)]);
Fiw(n)=norm([Qiw(1);Qfw(2)]);
Frw(n)=norm([Qrw(4);Qfw(5)]);


%Constraint error report
Phi=Phieval(tn,qn,par);
phiErr(n)=norm(Phi);
C=Ceval(tn,qn,par);
VelErr(n)=norm(C*qdn);
AccErr(n)=norm(C*qddn+Gam);

end

MaxKE=max(KE);
MinKE=min(KE);



   