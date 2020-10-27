%Tangent Space DAE Exp; Top Tip Fixed-Cardona Transient with gravity in
                        negative y-direction

utol=0.0000000001;     %Tolerance in solving for u
Btol=0.0000000001;     %Convergence criteria in B iteration
intol=0.000001;     %Tolerance in solving discretized equations of motion
Maxv=0.7;       %Limit on magnitude of v
Maxiter=6;      %Limit on number of integration iterations

h=0.0001;
h2=h^2;
tfinal=2;

g=9.81;
Iyy=0.234375;
Ixx=0.234375;
Izz=0.46875;
m=15;
par=[m;g;Ixx;Iyy;Izz];

integ=6;     %integ=1, Trap QuasiNR;  integ=2, SDIRK54 QuasiNR; 
             %integ=3, Nystrom4; 
             %integ=4, Trap PartNR; %integ=5, SDIRK54 PartNR
             %integ=6, Trap FullNR; %integ=7, SDIRK54 FullNR

% Data Storage Arrays
Q=zeros(7,10);
Qd=zeros(7,10);
Qdd=zeros(7,10);
Vv=zeros(3,10);
Vvd=zeros(3,10);
Vvdd=zeros(3,10);
LLam=zeros(4,10);
Uu=zeros(4,10);
Uud=zeros(4,10);
Uudd=zeros(4,10);
r=zeros(3,10);
Omegpr=zeros(3,10);
Rrpt=zeros(5,10);
Kk1=zeros(3,10);
Kk2=zeros(3,10);
Kk3=zeros(3,10);
Kk4=zeros(3,10);
Kk5=zeros(3,10);
R1=zeros(7,10);
R2=zeros(7,10);
R3=zeros(7,10);
R4=zeros(7,10);
R5=zeros(7,10);

% Initial Conditions

q0=[0;0;1;1;0;0;0];
p0=[q0(4);q0(5);q0(6);q0(7)];
kpr=[0;0;1];
omegapr0=[0;-4.61538;150];
pd0=0.5*GEval(p0)'*omegapr0;
rd0=atil(omegapr0)*kpr;
qd0=[rd0;pd0];
Q(:,1)=q0;
Qd(:,1)=qd0;

% Initial acceleration and Lagrange multiplier
QA=QAEval(q0,qd0,par);
S=SEval(q0,qd0,par);
[Gam,Gamsq,Gamsqd] = GamEval(q0,qd0);
M=MEval(q0,par);
Phiq = PhiqEval(q0);
Coef=[M,Phiq';Phiq,zeros(4,4)];
RHS=[S+QA;-Gam];
x=Coef\RHS;
qdd0=[x(1);x(2);x(3);x(4);x(5);x(6);x(7)];
Lam0=[x(8);x(9);x(10);x(11)];
Qdd(:,1)=qdd0;
LLam(:,1)=Lam0;

n=1;
npar=1;
t(1)=0;
vnormrpt(1)=Maxv+1;
jRepar=1;       %Counter for Reparameterization

jiterrpt(1)=Maxiter+1; %Maximum number of iterations in integration stages

while t(n)<tfinal;

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

vnorm=vnormrpt(n-1);
jiter=jiterrpt(n-1);

if vnorm>Maxv
Cr=2;       %Cr>1 triggers reparameterization

end

if jiter>Maxiter
Cr=Cr+2;
end

if Cr>1

npar=n-1;       %npar is last reparameterization step counter    
Crrpt(n)=Cr;
    
% Parameterization

[v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,Qdd,par,jRepar);

u=[0;0;0;0];
Uv(:,n-1)=u;
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;

jReparrpt(n)=jRepar;

%Jacobian Calculation

Lam0=LLam(:,n-1);
qd0=Qd(:,n-1);

[Rv,Rvd,Rvdd]=JacobEval(q0,qd0,vdd,Lam0,U,V,B,par);


Crr=1;
end

Cr=0;

% Integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if integ==1     %Tangent Space Trapezoidal QuasiNR %%%%%%%%%%%%%%%%%%%%%%%
    
    if Crr==1
Jvdd=Rvdd+(h/2)*Rvd+((h^2)/4)*Rv;
q=Q(:,n-1);
J=[Jvdd,PhiqEval(q)'];
Jcon(n)=cond(J);
Jinv=inv(J);
  
    Crr=0;
    end

[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapQuasiNR(n,npar,Vv,Vvd,Vvdd,...
    LLam,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par,Jinv);

R1nrpt(n)=R1n;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if integ==2     %Tangent Space SDIRK QuasiNR
    
    if Crr==1
% Integration Jacobian Evaluation
q=Q(:,n-1);
J=[Rvdd+h*0.25*Rvd+h2*0.0625*Rv,PhiqEval(q)'];
Jinv=inv(J);
Crr=0;

    end

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54QuasiNR(n,...
   npar,tn,Vv,Vvd,Vvdd,Qdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
   h,h2,utol,intol,Btol,par,Jinv);

R1(:,n)=RR(:,1);
R2(:,n)=RR(:,2);
R3(:,n)=RR(:,3);
R4(:,n)=RR(:,4);
R5(:,n)=RR(:,5);
Rrpt(:,n)=R;
Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
jiterrpt(n)=jiter;
hrpt(n)=h;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if integ==3      %Explicit Nystrom4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[v,vd,vdd,Lam,CondE,jiter]=IntegrateNystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,h2,utol,Btol,par);

Econdrpt(n)=CondE;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if integ==4     %Tangent Space Trapezoidal PartNR %%%%%%%%%%%%%%%%%%%%%%%%
    
    if Crr==1
Jvd=(h/2)*Rvd+((h^2)/4)*Rv;
  
    Crr=0;
    end

[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapPartNR(n,npar,Vv,Vvd,Vvdd,...
    LLam,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par,Jvd);

R1nrpt(n)=R1n;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if integ==5     %Tangent Space SDIRK PartNR
    
    if Crr==1
% Integration Jacobian Evaluation

Jvd=h*0.25*Rvd+h2*0.0625*Rv;

Crr=0;

    end

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54PartNR(n,...
   npar,tn,Vv,Vvd,Vvdd,Qdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
   h,h2,utol,intol,Btol,par,Jvd);

R1(:,n)=RR(:,1);
R2(:,n)=RR(:,2);
R3(:,n)=RR(:,3);
R4(:,n)=RR(:,4);
R5(:,n)=RR(:,5);
Rrpt(:,n)=R;
Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
jiterrpt(n)=jiter;
hrpt(n)=h;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if integ==6     %Tangent Space Trapezoidal FullNR %%%%%%%%%%%%%%%%%%%%%%%%

[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapFullNR(n,npar,Vv,Vvd,Vvdd,...
    Q,Qd,LLam,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par);

R1nrpt(n)=R1n;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if integ==7     %Tangent Space SDIRK FullNR
    

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54FullNR...
    (n,npar,tn,Vv,Vvd,Vvdd,Q,Qd,Qdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
    h,h2,utol,intol,Btol,par);

R1(:,n)=RR(:,1);
R2(:,n)=RR(:,2);
R3(:,n)=RR(:,3);    
R4(:,n)=RR(:,4);
R5(:,n)=RR(:,5);
Rrpt(:,n)=R;
Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
jiterrpt(n)=jiter;
hrpt(n)=h;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Process Results
Vv(:,n)=v;
Vvd(:,n)=vd;
Vvdd(:,n)=vdd;
LLam(:,n)=Lam;
vnormrpt(n)=norm(v);
vdnormrpt(n)=norm(vd);

jiterrpt(n)=jiter;

%Evaluate q

ue=Uu(:,n-1);
[u,Iteru] = usolv(ue,v,q0,V,U,B,utol);
Iterurpt(n)=Iteru;
Uu(:,n)=u;
unorm(n)=norm(u);
q=q0+V*v-U*u;
Q(:,n)=q;
qnormrpt(n)=norm(q);

%Update B
[B,Biter]=CorrectB(q,B,U,Btol);
Biterrpt(n)=Biter;

% Evaluate qd
Phiq = PhiqEval(q);
D=(eye(7)-U*B*Phiq)*V;
qd=D*vd;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

%Calculate Total Energy
r=[q(1);q(2);q(3)];
p=[q(4);q(5);q(6);q(7)];
kpr=[0;0;1];
M=MEval(q,par);
TEcalc=0.5*qd'*M*qd+par(1)*par(2)*q(2);
TE(n)=TEcalc;

% Evaluate qdd

P2=P2Eval(q,qd);
qdd=D*vdd-U*B*P2*qd;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

%Calculate quantities of interest
r(:,n)=[q(1);q(2);q(3)];
rx(n)=q(1);
ry(n)=q(2);
rz(n)=q(3);
p=[q(4);q(5);q(6);q(7)];
pd=[qd(4);qd(5);qd(6);qd(7)];
Omegpr(:,n)=2*GEval(p)*pd;
lam1(n)=Lam(1);
lam2(n)=Lam(2);
lam3(n)=Lam(3);
c=[15.234357;15.234357;0.46875];
Jpr=diag(c);
angmom0(:,n)=ATran(p)*Jpr*Omegpr(:,n);
uy=[0;1;0];
yandmonm0(n)=uy'*angmom0(:,n);

%Calculate constraint error
PosConstrNorm(n) = norm(PhiEval(q));
VelConstrNorm(n)=norm(Phiq*qd);
AccelErr=Phiq*qdd+P2*qd;
AccConstrNorm(n)=norm(AccelErr);

end


   