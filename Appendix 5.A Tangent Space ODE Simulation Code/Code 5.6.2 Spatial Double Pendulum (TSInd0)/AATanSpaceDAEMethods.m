%Tangent Space RK DAE; Spatial Double Pendulum

utol=0.000000001;     %Tolerance in solving for u
Btol=0.000000001;     %Convergence criteria in B iteration
intol=0.0001;     %Tolerance in solving discretized equations of motion
Maxv=1;       %Limit on magnitude of v
Maxiter=6;      %Limit on number of integration iterations
CondEmax=500;   %Limit on cond of coef matrix in Nystrom

h=0.001;
h2=h^2;
tfinal=20;

omeg1=0;
omeg2=5;

g=9.6;
Ixx=30;      %Ixx=Iyy=Izz
m=75;
par=[g;Ixx;m];

Integ=6;     %Integ=1, Trap QuasiNR;  Integ=2, SDIRK54 QuasiNR; 
            %Integ=3, Nystrom4; 
            %Integ=4, Trap PartNR; Integ=5, SDIRK54 PartNR
            %Integ=6, Trap FullNR; Integ=7, SDIRK54 FullNR

% Data Storage Arrays
Q=zeros(11,10);
Qd=zeros(11,10);
Qdd=zeros(11,10);
Vv=zeros(8,10);
Vvd=zeros(8,10);
Vvdd=zeros(8,10);
LLam=zeros(3,10);
Uu=zeros(3,10);
Uud=zeros(3,10);
Uudd=zeros(3,10);
Rrpt=zeros(5,10);
Kk1=zeros(8,10);
Kk2=zeros(8,10);
Kk3=zeros(8,10);
Kk4=zeros(8,10);
Kk5=zeros(8,10);

% Initial Conditions
q0=[1;0;0;0;1;0;0;0;0;0;-3];
p10=[q0(1);q0(2);q0(3);q0(4)];
p20=[q0(5);q0(6);q0(7);q0(8)];
r20=[q0(9);q0(10);q0(11)];
xpr=[1;0;0];
ypr=[0;1;0];
p1d0=(omeg1/2)*GEval(p10)'*xpr;
p2d0=(omeg2/2)*GEval(p20)'*ypr;
r2d0=[0;0;0];
qd0=[p1d0;p2d0;r2d0];

Q(:,1)=q0;
Qd(:,1)=qd0;
qnorm(1)=norm(q0);
qdnorm(1)=norm(qd0);

% Initial values of qdd and Lam
M=MEval(q0,par);
Phiq = PhiqEval(q0);
[Gam,Gamsq,Gamsqd] = GamEval(q0,qd0);
QA=QAEval(q0,qd0,par);
S=SEval(q0,qd0,par);
Coef=[M,Phiq';Phiq,zeros(3,3)];
Rhs=[QA+S;-Gam];
x=Coef\Rhs;
qdd=[x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10);x(11)];
Lam=[x(12);x(13);x(14)];
Qdd(:,1)=qdd;
LLam(:,1)=Lam;

n=1;
t(1)=0;
vnormrpt(1)=Maxv+1;
jRepar=0;       %Counter for Reparameterization
Maxjiterrpt(1)=Maxiter+1; %Maximum number of iterations in integration stages
jiterrpt(1)=0;

while t(n)<tfinal;

n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

vnorm=vnormrpt(n-1);
jiter=jiterrpt(n-1);

if vnorm>Maxv
Cr=2;

end

if jiter>Maxiter
Cr=Cr+2;
end

if Cr>1

npar=n-1;    
Crrpt(n)=Cr;
    
% Parameterization

[v,vd,vdd,q0,U,V,B,jRepar]=Param(n,Q,Qd,par,jRepar);

u=[0;0;0];
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;
jReparrpt(n)=jRepar;


%Jacobian Calculation

Lam=LLam(:,n-1);
qd0=Qd(:,n-1);
[Rv,Rvd,Rvdd]=JacobEval(q0,qd0,vdd,Lam,U,V,B,par);

Crr=1;
end

Cr=0;

% Integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==1     %Tangent Space Trapezoidal QuasiNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
   npar,tn,Vv,Vvd,Vvdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
   h,h2,utol,intol,Btol,par,Jinv);

Rrpt(:,n)=R;
Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
jiterrpt(n)=jiter;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==3      %Explicit Nystrom4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[v,vd,vdd,CondE,jiter]=IntegrateNystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,h2,utol,Btol,par);

Econdrpt(n)=CondE;

if CondE>CondEmax
    Cr=2;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==4     %Tangent Space Trapezoidal PartialNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
% Integration Jacobian Evaluation

Jvd=h*0.25*Rvd+h2*0.0625*Rv;

Crr=0;

    end

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54PartNR(n,...
   npar,tn,Vv,Vvd,Vvdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
   h,h2,utol,intol,Btol,par,Jvd);

Rrpt(:,n)=R;
Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
jiterrpt(n)=jiter;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==6     %Tangent Space Trapezoidal PartialNR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[v,vd,vdd,Lam,jiter,R1n]=IntegTanSpTrapFullNR(n,npar,Vv,Vvd,Vvdd,...
    LLam,Uu,q0,V,U,B,h,h2,utol,intol,Btol,par);

R1nrpt(n)=R1n;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Integ==7     %Tangent Space SDIRK FullNR

[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54FullNR(n,...
   npar,tn,Vv,Vvd,Vvdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
   h,h2,utol,intol,Btol,par);

Rrpt(:,n)=R;
Kk1(:,n)=k1;
Kk2(:,n)=k2;
Kk3(:,n)=k3;
Kk4(:,n)=k4;
Kk5(:,n)=k5;
jiterrpt(n)=jiter;

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
Phiq = PhiqEval(q);
[B,Biter]=CorrectB(q,B,U,Btol);
Biterrpt(n)=Biter;

% Evaluate qd

D=(eye(11)-U*B*Phiq)*V;
qd=D*vd;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);

%Calculate Total Energy
kpr=[0;0;1];
M=MEval(q,par);
TEcalc=0.5*qd'*M*qd+par(1)*par(3)*q(11);
TE(n)=TEcalc;

% Evaluate qdd

P2=P2Eval(q,qd);
qdd=D*vdd-U*B*P2*qd;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);

%Calculate quantities of interest
r2x(n)=q(9);
r2y(n)=q(10);
r2z(n)=q(11);

%Calculate constraint error
PosConstrNorm(n) = norm(PhiEval(q));
VelConstrNorm(n)=norm(Phiq*qd);
AccelErr=Phiq*qdd+P2*qd;
AccConstrNorm(n)=norm(AccelErr);


end


   