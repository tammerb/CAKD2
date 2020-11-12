
%General Tangent Space ODE-Planar Double Pendulum Example

utol=10^-8;     %Tolerance in solving for u
Btol=10^-8;     %Convergence criteria in B iteration
intol=10^-6;    %Tolerance in solving discretized equations of motion
Atol=10^-6;     %Absolute error tolerance for variable step metnods
Maxv=1;         %Limit on magnitude of v
MaxEcond=300;     %Limit on magnitude of reduced mass matrix condition number
Maxiter=8;      %Limit on number of implicit integration iterations

h=0.001;        %Step size
hmax=0.001;     %Maximum Allowable Step size
hvar=1;         %hvar=1, variable step;hvar=2, fixed step

tfinal=9.999;     %Final time

Integ=1;    %Numerical Integration Options

    %Explicit Integration Methods:
        %Integ=1, Nystrom4; Integ=2, RungeKutta4;
        %Integ=3, Kutta3/8; integ=4, RKFN
                                    
    %Implicit Integration Methods:
        %Integ=5, Trapezoidal; Integ=6, SDIRK54
        
% Variable dimension data
nq=6;
nh=4;

% Inertia and Force Data
m1=1;
m2=1;
Jp1=0.3;
Jp2=0.3;
g=9.8;
K1=20;
K2=20;
K3=0;    %There is no K3 spring   
C1=10;
C2=10;
n1=0;
n2=0;

%Fixed Parameter data List-Complete Partitioning in function parPart
par=[nq;nh;utol;Btol;intol;Atol;m1;m2;Jp1;Jp2;g;K1;K2;K3;C1;C2;n1;n2];

% Data Storage Arrays
Q=zeros(nq,10);
Qd=zeros(nq,10);
Qdd=zeros(nq,10);
Vv=zeros(nq-nh,10);
Vvd=zeros(nq-nh,10);
Vvdd=zeros(nq-nh,10);
Uu=zeros(nh,10);

% Initial Conditions
q0=[0;-1;-pi/2;0;-3;-pi/2];
qd0=[0;0;0;10;0;10];
Q(:,1)=q0;
phi1(1)=q0(3);
phi2(1)=q0(6);
Qd(:,1)=qd0;

%Integration Preparation
n=1;
t(1)=0;

%Criteria to Assure Initial Parameterization
normv=Maxv+1;
Iterinteg=0;
irepr=0;
TEcalc=0;

jRepar=1;
Econdrpt(1)=MaxEcond+1;
vnorm(1)=10;
jiterrpt(1)=Maxiter+1;
nch=1;
rpth(1)=h;

% Integration

while t(n)<tfinal
%Time Step Update
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

%Parameterization
if Integ<5  %Reparameterization Criteria for Explicit Integration 
    
if vnorm(n-1)>Maxv
    Cr=2;
end
if Econdrpt(n-1)>MaxEcond
    Cr=Cr+2;
end
end

if Integ>4  %Reparameterization Criteria for Implicit Integration
    
if vnorm(n-1)>Maxv
Cr=2;
end
if jiterrpt(n-1)>Maxiter
Cr=Cr+2;
end          
end

if Cr>1     %Criteria for reparameterization
 
Crrpt(n)=Cr;

% Parameterization
[vd,q0,U,V,B,jRepar]=Param(t,n,Q,Qd,jRepar,par);

u=zeros(nh,1);
v=zeros(nq-nh,1);
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Uu(:,n-1)=u;
jReparrpt(n)=jRepar;

Cr=0;
end

% Integration
tnm=t(n-1);
qnm=Q(:,n-1);
qdnm=Qd(:,n-1);

vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
unm=Uu(:,n-1);
ue=unm;

if Integ<5      %Explicit integrators
if Integ==1
[vn,vdn,vddn,Econd]=ExplicitNystrom4(tnm,vnm,vdnm,ue,...
    U,V,B,q0,par,h);

end

if Integ==2
    [vn,vdn,vddn,Econd]=ExplicitRungeKutta4(tnm,vnm,vdnm,ue,...
    U,V,B,q0,par,h);
end

if Integ==3
    [vn,vdn,vddn,Econd]=ExplicitKutta38(tnm,vnm,vdnm,ue,...
    U,V,B,q0,par,h);
end

if Integ==4
[vn,vdn,vddn,Econd,h,nch]=ExplicitRKFN45(n,tnm,vnm,vdnm,ue,...
    U,V,B,q0,h,hmax,par,Atol,nq,nh,nch,hvar);
hrpt(n)=h;
end

Econdrpt(n)=Econd;
end

if Integ >4     %Implicit integrators
    
if Integ==5
[vn,vdn,vddn,jiter,R1Norm,JCond]=ImplicitTrap(n,tn,...
    Vv,Vvd,Vvdd,Uu,Q,Qd,Qdd,q0,U,V,B,intol,par,h);    
        
end
if Integ==6
    
[vn,vdn,vddn,jiter,R1Norm,JCond,h,err]=ImplicitSDIRK54(n,tn,...
    Vv,Vvd,Vvdd,Uu,Q,Qd,Qdd,q0,U,V,B,par,intol,Atol,nq,nh,h,hmax,nch,hvar); 
hrpt(n)=h;
R1Normrpt(n)=R1Norm;
errrpt(n)=err;
    
end

JCondrpt(n)=JCond;
jiterrpt(n)=jiter;

end

%Evaluate and Record Solution
Vv(:,n)=vn;
vnorm(n)=norm(vn);
Vvd(:,n)=vdn;
vdnorm(n)=norm(vdn);
Vvdd(:,n)=vddn;
vddnorm(n)=norm(vddn);

%Evaluate/update un, qn and B
[u,Iteru]=ueval(tn,ue,vn,q0,V,U,B,par);
un=u;
Uu(:,n)=un;
unorm(n)=norm(u);
uIterrpt(n)=Iteru;
qn=q0+V*vn-U*un;
Q(:,n)=qn;
qnorm(n)=norm(qn);

[B,Biter]=Beval(tn,qn,B,U,par);
Biterrpt(n)=Biter;

%Evaluate qdn and qddn

I=eye(nq);
D=Deval(tn,qn,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,qn,par);
qdn=D*vdn-U*B*Pst;
Qd(:,n)=qdn;
qdnorm(n)=sqrt(qdn'*qdn);
qddn=D*vddn-U*B*Gameval(tn,qn,qdn,par);
Qdd(:,n)=qddn;
qddnorm(n)=norm(qddn);

%Report key data
trpt(n)=tn;
x1(n)=qn(1);
y1(n)=qn(2);
phi1(n)=qn(3);
x2(n)=qn(4);
y2(n)=qn(5);
phi1(n)=qn(3);
phi2(n)=qn(6);
phi2m1(n)=qn(6)-qn(3);
phi2m1d(n)=qdn(6)-qdn(3);
phi2m1dd(n)=qddn(6)-qddn(3);

%Calculate Total Energy

M=Meval(qn,par);
TEcalc=0.5*(qdn'*M*qdn)+g*(qn(2)+qn(5))+...
    0.5*(K1*(qn(3)+pi/2)^2+K2*(qn(6)-qn(3))^2);
TE(n)=TEcalc;


%Constraint error calculation
q=Q(:,n);
qd=Qd(:,n);
qdd=Qdd(:,n);
PosConstrNorm(n) = norm(Phi(tn,qn,par));
Phiq = P1(t,q,par);
VelConstrNorm(n)=norm(Phiq*qd);
AccelErr=Phiq*qdd+P2(tn,q,qd,par)*qd;
AccConstrNorm(n)=norm(AccelErr);



end




   