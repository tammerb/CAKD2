%AA Planar Tangent Space Multibody Simulation

%Integration and Error Control Parameters
utol=10^-12;         %Tolerance in solving for u
Btol=10^-12;         %Convergence criteria in B iteration
intol=10^-6;        %Tolerance in solving discretized equations of motion
Atol=10^-6;         %Absolute error tolerance for variable step metnods
Maxv=10;            %Limit on norm of v
Maxu=10;            %Limit on norm of u
MaxImpSoliter=10;    %Limit on number of implicit integration iterations
MaxUiter=8;         %Limit on number of U iterations
MaxJcond=5000;       %Limit on magnitude of Jcond in Trap
R1nmax=30000;        %Limit on initial residual in implicit integration
MaxECond=200;        %Limit on magnitude of ECond
MaxBCond=600000;      %Limit on magnitude of BCond
MaxBnormRat=10;     %Maximum ratio of Bnorm to Bnorm at parameterization
MaxBCondRat=10;     %Maximum ratio of BCond to BCond0 at parameterization

h0=0.001;           %Initial time step
h=h0;
hmax=0.001;          
hvar=1;             %hvar=1, variable h;hvar=2, constant h

tfinal=2;           %Final simulation time

integ=1;            %integ=1-Impl ODE Trap; integ=2-Impl ODE SDIRK54; 
                    %integ=3-Impl Ind0 Trap; integ=4-Impl Ind0 SDIRK54;
                    %integ=5-Expl ODE Nystrom4; integ=6-Expl ODE RKFN45;
                    %integ=7-Expl Ind0 Nystrom4; integ=8-Expl Ind0 RKFN45
                    
InvJ=1;             % InvJ=1, Invert J; InvJ=2, Factor J
            
g=9.8;

%Application Data
    %app=1, Double Pendulum
    %app=2, Quick Return
    %app=3, Lumped Mass Coil Spring-5 masses
    %app=4, Lumped Mass Coil Spring-10 masses
    %app=5, Three Body Translational Model
    %app=6, Slider-Crank
    %app=7, Rotating Disk with Translating Body
    %app=8, Multiple Slider-Crank
    %app=9, Flywheel-Spring
    %app=10, Loader
    %app=11, Slider-Crank-2
    
app=1;

[nb,ngc,nh,nc,nv,nu,NTSDA,NRSDA,PJDT,PMDT,PTSDAT,PRSDAT,q0,qd0]=...
    AppData(app);     %AppData(app)

par=[nb;ngc;nh;nc;nv;nu;g;utol;Btol;intol;Atol;h0;hvar;NTSDA;NRSDA];

%Initial condition calculation, if required
if app==2       
phi2d0=6;
Phiq=PhiqEval(0,q0,PJDT,par);
EEE=[Phiq;zeros(1,5),1,zeros(1,6)];
RRHS=[zeros(11,1);phi2d0];
qd0=EEE\RRHS;   
end

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);
 
% Data Storage Arrays
Q=zeros(ngc,10);
Qd=zeros(ngc,10);
Qdd=zeros(ngc,10);
Vv=zeros(nv,10);
Vvd=zeros(nv,10);
Vvdd=zeros(nv,10);
Uu=zeros(nu,10);
Uud=zeros(nu,10);
Uudd=zeros(nu,10);
LLam=zeros(nc,10);

Q(:,1)=q0;
Qd(:,1)=qd0;

%Calculate Initial Acceleration and Lagrange Multipliers
Phiq=PhiqEval(0,q0,PJDT,par);
Phiq0=Phiq;
U=Phiq';
B=inv(U'*U);
V=null(U');
v0=zeros(nv,1);
vd0=V'*qd0;
u0=zeros(nu,1);
Vv(:,1)=v0;
Vvd(:,1)=vd0;
Uu(:,1)=u0;

[vdd0,ECond]=ODEfunct(0,v0,vd0,PMDT,PTSDAT,PRSDAT,PJDT,...
    u0,q0,V,U,B,par);
D=(eye(ngc)-U*B*Phiq)*V;
Gam=GamEval(0,q0,qd0,PJDT,par);
qdd0=D*vdd0-U*B*Gam;
Lam0=B'*U'*(MEval(PMDT,par)*qdd0+QAEval(0,q0,qd0,PMDT,PTSDAT,PRSDAT,par));

Qdd(:,1)=qdd0;
LLam(:,1)=Lam0;
Vvdd(:,1)=vdd0;

%Initialize Data For Integration
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
jRepar=0;       %Counter for Reparameterization time step
Cr=2;           %Reparameterize if Cr>1
jiter=0;
Iterurpt(1)=0;
JCondrpt(1)=0;
ImpSoliterrpt(1)=0;
R1Normrpt(1)=0;
vnormrpt(1)=0;
unormrpt(1)=0;
ECondrpt(1)=0;
Econd(1)=1;
uiterrpt(1)=1;
nch=1;
VelConstrNorm(1)=0;
AccConstrNorm(1)=0;
nres=1;
jnres=0;
B0norm=1;
Bnorm(1)=1;
BCond(1)=MaxBCond+1;
BCondRat(1)=1;
Jinv=eye(nv);

%Integration
while t(n)<tfinal
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

%Reparameterization Criteria
if integ<5     %Start implicit reparameterization criteria 
if BCondRat(n-1)>10
Cr=2; 
end
if vnormrpt(n-1)>Maxv
Cr=2;
end
if unormrpt(n-1)>Maxu
Cr=2;
end
if ImpSoliterrpt(n-1)>MaxImpSoliter
Cr=Cr+2;
end
if JCondrpt(n-1)>MaxJcond
Cr=Cr+2;
end
if R1Normrpt(n-1)>R1nmax
    Cr=Cr+2;
end
if Bnorm(n-1)/B0norm>MaxBnormRat
    Cr=Cr+2;
end  
if BCond(n-1)>MaxBCond
    Cr=Cr+2;
end
end      %end implicit reparameterization criteria

if integ>4      %Start explicit reparameterization criteria    
if vnormrpt(n-1)>Maxv
Cr=2;
end
if unormrpt(n-1)>Maxu
Cr=2;
end    
if ECondrpt(n-1)>MaxECond
Cr=Cr+2;
end
if Bnorm(n-1)/B0norm>MaxBnormRat
    Cr=Cr+2;
end
if BCond(n-1)>MaxBCond
    Cr=Cr+2;
end    
end            %end explicit reparameterization criteria

if Cr>1     %Reparameterization
npar=n-1;    
Crrpt(n)=Cr;
nch=1;
tnm=t(n-1);
[vd,vdd,q0,U,V,B,jRepar]=Param(n,tnm,Q,Qd,Qdd,PJDT,par,jRepar);
u=zeros(nu,1);
v=zeros(nv,1);
Uu(:,n-1)=u;
Vv(:,n-1)=v;
Vvd(:,n-1)=vd;
Vvdd(:,n-1)=vdd;
B0norm=norm(B);
BCond0=B0norm*norm(U'*U);
Cr=1;
jReparrpt(n)=jRepar;
end        %End Reparameterization

% Integration
if integ==1     %Implicit ODE trapezoidal    
[v,vd,vdd,R1n,ImpSoliter,JCond,Jinv,Jinviter,h,nch]=...
    ImplicitODETrap(n,tn,npar,Vv,Vvd,Vvdd,Uu,q0,V,U,B,...
    h,hmax,nch,PMDT,PTSDAT,PRSDAT,PJDT,par,InvJ,Jinv);
R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
ImpSoliterrpt(n)=ImpSoliter;
JinvIterrpt(n)=Jinviter;
hrpt(n)=h;
end

if integ==2     %Implicit ODE SDIRK54 
[v,vd,vdd,ImpSoliter,R1n,J,JCond,h,nch,Jinv,Jinviter]=...
    ImplicitODESDIRK54(n,tn,Vv,Vvd,Vvdd,Uu,Q,Qd,Qdd,q0,U,V,B,par,...
    h,hmax,nch,npar,PJDT,PMDT,PTSDAT,PRSDAT,Jinv,InvJ);
R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
ImpSoliterrpt(n)=ImpSoliter;
hrpt(n)=h;
Jinviterrpt(n)=Jinviter;
end

if integ==3        %Implicit Index 0 trapezoidal
[v,vd,vdd,Lam,R1n,ImpSoliter,JCond,h,nch]=...
    ImplicitInd0Trap(n,tn,Vv,Vvd,Vvdd,Uu,LLam,q0,V,U,B,...
    h,hmax,nch,PMDT,PTSDAT,PRSDAT,PJDT,par);
R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
ImpSoliterrpt(n)=ImpSoliter;
hrpt(n)=h;    
end

if integ==4     %Implicit Index 0 SDIRK54
[v,vd,vdd,Lam,ImpSoliter,R1n,J,JCond,h,nch]=...
    ImplicitInd0SDIRK54(n,tn,Vv,Vvd,Vvdd,LLam,Uu,q0,U,V,B,par,...
    h,hmax,nch,npar,PJDT,PMDT,PTSDAT,PRSDAT);
R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
ImpSoliterrpt(n)=ImpSoliter;
hrpt(n)=h;
end

if integ==5     %Explicit ODE Nystrom 
[v,vd,vdd,ECond]=ExplicitODENystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,npar,PMDT,PTSDAT,PRSDAT,PJDT,par);
ECondrpt(n)=ECond;    
end

if integ==6     %Explicit ODE RKFN45  
[v,vd,vdd,ECond,Err,h,nch]=ExplicitODERKFN45(n,tn,Vv,Vvd,Uu,...
    U,V,B,q0,h,hmax,par,PMDT,PTSDAT,PRSDAT,PJDT,nch);
ECondrpt(n)=ECond;
Errrpt(n)=Err;
hrpt(n)=h;
end

if integ==7         %Explicit Index0 Nystrom4
[v,vd,vdd,Lam,ECond]=ExplicitInd0Nystrom4(n,tn,Vv,Vvd,Uu,V,U,B,q0,...
    h,npar,PMDT,PTSDAT,PRSDAT,PJDT,par); 
ECondrpt(n)=ECond;
LLam(:,n)=Lam;
end

if integ==8         %Explicit Index0 RKFN45
[v,vd,vdd,Lam,ECond,Err,h,nch]=...
    ExplicitInd0RKFN45(n,tn,Vv,Vvd,Uu,...
    U,V,B,q0,h,hmax,par,PMDT,PTSDAT,PRSDAT,PJDT,nch);
ECondrpt(n)=ECond;
Errrpt(n)=Err;
hrpt(n)=h;
LLam(:,n)=Lam;
end

%Process Results
Vv(:,n)=v;
Vvd(:,n)=vd;
Vvdd(:,n)=vdd;
vnormrpt(n)=norm(v);
vdnormrpt(n)=norm(vd);
vddnormrpt(n)=norm(vdd);

%Evaluate q
u=Uu(:,n-1);
if n-1>npar
u=2*u-Uu(:,n-2);
end
[u,Iteru]=usolv(tn,u,v,q0,PJDT,V,U,B,par);
Iterurpt(n)=Iteru;
Uu(:,n)=u;
unormrpt(n)=norm(u);
q=q0+V*v-U*u;
Q(:,n)=q;
qnormrpt(n)=norm(q);

%Evaluate B and qd
[B,Biter]=BEval(tn,q,B,U,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
D=(eye(ngc)-U*B*Phiq)*V;
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
qd=D*vd-U*B*Pst;
Qd(:,n)=qd;
qdnormrpt(n)=norm(qd);
Bnorm(n)=norm(B);
BnormRat(n)=Bnorm(n)/B0norm;
Biterrpt(n)=Biter;
BCond(n)=Bnorm(n)*norm(Phiq*U);
BCondRat(n)=BCond(n)/BCond0;

% Evaluate qdd and Lam (Lam evaluation a postprocessing step)
Gam=GamEval(tn,q,qd,PJDT,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;
qddnormrpt(n)=norm(qdd);
Lam=B'*U'*(-MEval(PMDT,par)*qdd+QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par));
LLam(:,n)=Lam;

%End integration

%Calculate output data

%Total Energy
%Kinetic Energy
M=MEval(PMDT,par);
KE(n)=0.5*qd'*M*qd;

%Potential Energy of Gravitational Force
PEg=0;
i=1;
while i<=nb
m=PMDT(1,i);
PEg=PEg+m*g*q(2+3*(i-1));
i=i+1;
end

%Potential Energy of Spring Force
PEs=0;
if NTSDA>0
T=1;
while T<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=PTSDATPart(PTSDAT,T);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
el=sqrt(d12'*d12);

%K=K*(1-sign(q(1)))/2;     %Unilateral spring constant
PEs=PEs+0.5*K*(el-el0)^2;
T=T+1;
end

end

if NRSDA>0
R=1;
while R<=NRSDA
[i,j,K,C,phi0,T]=PRSDATPart(PRSDAT,R);
[r1,ph1]=qPart(q,i);
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
PEs=PEs+0.5*K*(ph2-ph1-phi0)^2;
R=R+1;
end
end

PE(n)=PEg+PEs;

TE(n)=KE(n)+PE(n);

%Data of Interest (Enter for each application)
if app==1   %Double Pendulum
x1(n)=q(1); 
y1(n)=q(2); 
phi1(n)=q(3); 
x2(n)=q(4); 
y2(n)=q(5); 
phi2(n)=q(6);
phi2m1(n)=phi2(n)-phi1(n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Quick Return
phi1(n)=q(3);
phi2(n)=q(6);
x4(n)=q(10);
x4d(n)=qd(10);
x4dd(n)=qdd(10);
Lam6(n)=Lam(6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Lumped Mass Coil Spring-5 Masses 
x1del(n)=q(1)-0.2; 
x2del(n)=q(4)-0.4;
x3del(n)=q(7)-0.6;
x4del(n)=q(10)-0.8;
x5del(n)=q(13)-1;

x1d(n)=qd(1); 
x2d(n)=qd(4);
x3d(n)=qd(7);
x4d(n)=qd(10);
x5d(n)=qd(13);

x1dd(n)=qd(1); 
x2dd(n)=qd(4);
x3dd(n)=qd(7);
x4dd(n)=qd(10);
x5dd(n)=qd(13);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4  %Lumped Mass Coil Spring-10 Masses
x1del(n)=q(1)-0.1; 
x2del(n)=q(4)-0.2;
x3del(n)=q(7)-0.3;
x4del(n)=q(10)-0.4;
x5del(n)=q(13)-0.5;
x6del(n)=q(16)-0.6; 
x7del(n)=q(19)-0.7;
x8del(n)=q(22)-0.8;
x9del(n)=q(25)-0.9;
x10del(n)=q(28)-1;

x1d(n)=qd(1); 
x2d(n)=qd(4);
x3d(n)=qd(7);
x4d(n)=qd(10);
x5d(n)=qd(13);
x6d(n)=qd(16); 
x7d(n)=qd(19);
x8d(n)=qd(22);
x9d(n)=qd(25);
x10d(n)=qd(28);

x1dd(n)=qd(1); 
x2dd(n)=qd(4);
x3dd(n)=qd(7);
x4dd(n)=qd(10);
x5dd(n)=qd(13);
x6dd(n)=qd(16); 
x7dd(n)=qd(19);
x8dd(n)=qd(22);
x9dd(n)=qd(25);
x10dd(n)=qd(28);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5   %Three Body Translational Model
y1(n)=q(2); 
x2(n)=q(4); 
x3(n)=q(7);
y1d(n)=qd(2); 
x2d(n)=qd(4); 
x3d(n)=qd(7);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==6   %Slider-Crank
phi1(n)=q(3);
phi1d(n)=qd(3);
phi1dd(n)=qdd(3);
x2(n)=q(4);
x2d(n)=qd(4);
Lam=LLam(:,n);
Lam1=[Lam(1);Lam(2)];
FRev(n)=norm(Lam1)/10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if app==7  %Rotating Disk with Translating Body
phi1rpt(n)=q(3);
A1=ATran(q(3));
r2=[q(4);q(5)];
dely2pr(n)=(A1'*(r2-A1*ux))'*uy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if app==8   %Multiple Slider-Crank
phi1(n)=q(3);
phi1d(n)=qd(3);
x2(n)=q(4);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==9   %Flywheel-Spring
phi1(n)=q(3);
phi1d(n)=qd(3);
Lam1(n)=Lam(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==10   %Loader
y1(n)=q(2);
y2(n)=q(5);
phi1(n)=q(3);
phi2(n)=q(6);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate constraint error
Phi=PhiEval(tn,q,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
Gam=GamEval(tn,q,qd,PJDT,par);
PosConstrNorm(n)=norm(Phi);
VelConstrNorm(n)=norm(Phiq*qd);
AccConstrNorm(n)=norm(Phiq*qdd+Gam);


end

