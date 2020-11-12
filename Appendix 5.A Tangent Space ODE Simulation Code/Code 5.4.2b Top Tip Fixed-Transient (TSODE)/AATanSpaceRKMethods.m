
%5.4.2b Top Tip Fixed-Cardona Transient

%Functions the user must enter
    %Meval,M2
    %Phi, P1, P2, P3, P4, P5, P6
    %QAeval, QAsqqd, Seval, Ssqqd
%Functions the user may choose to enter
    %parPart
    %qPart, qdPart

utol=10^-12;    %Tolerance in solving for u
Btol=10^-12;    %Convergence criteria in B iteration
intol=10^-6;    %Tolerance in solving discretized equations of motion
Atol=10^-6;     %Absolute error tolerance for variable step metnods
Maxv=0.75;       %Limit on magnitude of v
MaxEcond=50;    %Limit on magnitude of reduced mass matrix condition number
Maxiter=8;      %Limit on number of implicit integration iterations
hvar=1;         %hvar=1, variable step;  hvar=2, fixed step
hmax=0.01;      %Maximum h in varible step size integraation

h=0.0001;

tfinal=2;

Integ=5;       %Integration method: Integ=1, Nystrom; Integ=2, RungeKutta4;
                                    %Integ=3, Kutta3/8  %Integ=4, RKFN45
                                    %Integ=5, Trap; Integ=6, SDIRK54
nq=7;
nh=4;
% Inertia and Force Data, Add as Needed
g=9.81;
m=15;
Ixy=0.234357;
Iz=0.46875;
Dcf=0;              %Aerodynamic drag coefficient

%Fixed Parameter data List-Complete Partitioning in function parPart
par=[nq;nh;utol;Btol;intol;Atol;m;g;Ixy;Iz;Dcf];

%Global Unit Vectors
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

% Data Storage Arrays
Q=zeros(nq,10);
Qd=zeros(nq,10);
Qdd=zeros(nq,10);
Vv=zeros(nq-nh,10);
Vvd=zeros(nq-nh,10);
Vvdd=zeros(nq-nh,10);
Uu=zeros(nh,10);
rRpt=zeros(3,10);
Omegapr=zeros(3,10);

% Initial Conditions
r0=[0;0;1];
p0=[1;0;0;0];
q0=[r0;p0];
Omegpr0=[0;-4.61538;150];
rd0=atileval(Omegpr0)*uz;
pd0=0.5*Geval(p0)'*Omegpr0;
qd0=[rd0;pd0];
Q(:,1)=q0;
Qd(:,1)=qd0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO USER CHANGES/INPUT REQUIRED BEYOND THIS POINT

%Integration Preparation
n=1;
t(1)=0;

normv=Maxv+1;
Iterinteg=0;
irepr=0;
TEcalc=0;

jRepar=1;
Econdrpt(1)=MaxEcond+1;
vnorm(1)=10;
jiterrpt(1)=Maxiter+1;
nch=1;

% Integration
while t(n)<tfinal
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

[vd,q0,U,V,B,jRepar]=Param(t,n,Q,Qd,jRepar,par);

u=zeros(nh,1);
v=zeros(nq-nh,1);
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Uu(:,n-1)=u;
jReparrpt(n)=jRepar;

Cr=0;
end

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
[r,p]=qPart(qn);
[rd,pd]=qdPart(qdn);
rRpt(:,n)=r;
rx(n)=r(1);
ry(n)=r(2);
ryd(n)=rd(2);
rz(n)=r(3);
Omegapr(:,n)=2*Geval(p)*pd;
OmegaprNorm(n)=norm(Omegapr(:,n));
Jpr=diag([Ixy;Ixy;Iz]);
zprAngMomNorm(n)=norm((Jpr*2*Geval(p)*pd)'*uz);

%Calculate Kinetic Energy; Add Potential Energy for Total Energy

M=Meval(qn,par);
KE(n)=0.5*qdn'*M*qdn;
TE(n)=KE(n)+m*g*ry(n);


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




   