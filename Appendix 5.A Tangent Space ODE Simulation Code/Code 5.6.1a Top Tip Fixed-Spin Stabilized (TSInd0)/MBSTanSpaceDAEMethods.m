%Tangent Space RK DAE Exp; Top Tip Fixed

utol=0.0000000001;     %Tolerance in solving for u
Btol=0.0000000001;     %Convergence criteria in B iteration
intol=0.000001;     %Tolerance in solving discretized equations of motion
Maxv=0.7;       %Limit on magnitude of v
Maxiter=6;      %Limit on number of integration iterations

h=0.001;
h2=h^2;
tfinal=1;

g=9.6;
Ixy=90;      %Ixx=Iyy=Ixy
Iz=30;
m=30;
omeg=12.5;
par=[m;g;Ixy;Iz];

Integ=4;  %Integ=1, Trap QuasiNR;  Integ=2, SDIRK54 QuasiNR; 
         %Integ=3, Nystrom4; 
         %Integ=4, Trap PartNR; Integ=5, SDIRK54 PartNR
         %Integ=6, Trap FullNR; Integ=7, SDIRK54 FullNR

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
mom=zeros(3,1);
angmom0=zeros(3,1);

% Initial Conditions

epsx=0.000000000001;        %Initial x- and y- angular velocity
epsy=0.000000000001;
q0=[0;0;1;1;0;0;0];
p0=[q0(4);q0(5);q0(6);q0(7)];
kpr=[0;0;1];
omegapr0=[epsx;epsy;omeg];
pd0=0.5*GEval(p0)'*omegapr0;
rd0=ATran(p0)*atil(omegapr0)*kpr;
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
jRepar=0;       %Counter for Reparameterization

jiterrpt(1)=Maxiter+1; %Maximum number of iterations in integration stages

while t(n)<tfinal;

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

jiter=jiterrpt(n-1);

if vnormrpt(n-1)>Maxv
Cr=2;

end

if jiter>Maxiter
Cr=Cr+2;
end

if Cr>1

npar=n-1;    
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

Lam=LLam(:,n-1);
qd=Qd(:,n-1);

[Rv,Rvd,Rvdd]=JacobParam(q0,qd,vdd,Lam,U,V,B,par);

%[ErrRvdd,ErrRvd,ErrRv,RvddEst,RvdEst,RvEst] = JacobianCheck(v,...
    %vd,vdd,Lam0,u,par,q0,V,U,B,utol,Btol,h,Rvdd,Rvd,Rv);

Crr=1;
end

Cr=0;

% Integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==1     %Tangent Space Trapezoidal QuasiNR %%%%%%%%%%%%%%%%%%%%%%%
    
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
if Integ==2     %Tangent Space SDIRK QuasiNR
    
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
if Integ==3      %Explicit Nystrom4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[v,vd,vdd,Lam,CondE,jiter]=IntegrateNystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,h2,utol,Btol,par);

Econdrpt(n)=CondE;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==4     %Tangent Space Trapezoidal PartNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Crr==1
Jvd=(h/2)*Rvd+((h^2)/4)*Rv;

  
    Crr=0;
    end

[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapPartNR(n,npar,Vv,Vvd,Vvdd,...
    LLam,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par,Jvd);

R1nrpt(n)=R1n;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==5     %Tangent Space SDIRK PartNR
    
    if Crr==1
% Integration SubJacobian Evaluation

Jvd=h*0.25*Rvd+h2*0.0625*Rv;

Crr=0;

    end

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54PartNR...
    (n,npar,tn,Vv,Vvd,Vvdd,Qdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
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
if Integ==6     %Tangent Space Trapezoidal FullNR %%%%%%%%%%%%%%%%%%%%%%%
    

[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapFullNR(n,npar,Vv,Vvd,Vvdd,...
    LLam,Q,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par);

R1nrpt(n)=R1n;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==7     %Tangent Space SDIRK FullNR

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54FullNR...
    (n,npar,tn,Vv,Vvd,Vvdd,Q,Qdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
TEcalc=0.5*qd'*M*qd+par(1)*par(2)*q(3);
TE(n)=TEcalc;

% Evaluate qdd

[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
qdd=D*vdd-U*B*Gam;
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
rd=[qd(1);qd(2);qd(3)];
mom(:,n)=m*rd;
c=[120;120;30];
Jpr=diag(c);
angmom0(:,n)=ATran(p)*Jpr*Omegpr(:,n);
uz=[0;0;1];
zangmom0(n)=uz'*angmom0(:,n);

%Calculate constraint error
PosConstrNorm(n) = norm(PhiEval(q));
VelConstrNorm(n)=norm(Phiq*qd);
AccelErr=Phiq*qdd+Gam;
AccConstrNorm(n)=norm(AccelErr);

end


   