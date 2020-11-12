%AA Planar Index 0 Multibody Simulation with Friction

%Integration and Error Control Parameters
utol=10^-12;         %Tolerance in solving for u
Btol=10^-12;         %Convergence criteria in B iteration
intol=10^-6;        %Tolerance in solving discretized equations of motion
Atol=10^-5;         %Position error tolerance for variable step metnods
Vtol=10^-4;         %Velocity error tolerance for variable step metnods
Maxv=1;            %Limit on norm of v
Maxu=10;            %Limit on norm of u
MaxImpSoliter=12;    %Limit on number of implicit integration iterations
MaxUiter=8;         %Limit on number of U iterations
MaxJcond=50000;       %Limit on magnitude of Jcond in Trap
R1nmax=3000;        %Limit on initial residual in implicit integration
MaxECond=50000;        %Limit on magnitude of ECond
MaxBCond=10^6;      %Limit on magnitude of BCond

MaxBnormRat=10;     %Maximum ratio of Bnorm to Bnorm at parameterization
MaxBCondRat=10;     %Maximum ratio of BCond to BCond0 at parameterization

h0=0.001;           %Initial time step
h=h0;
hmax=0.001;          
hvar=1;             %hvar=1, variable h;hvar=2, constant h

tfinal=4;           %Final simulation time

integ=1;            %integ=1-Impl Ind0 Trap; integ=2-Impl Ind0 SDIRK54;
                    %integ=3-Expl Ind0 Nystrom4; integ=4-Expl Ind0 RKFN45;
                    %integ=5-Expl Ind0 Forward Euler
            
g=0;              %gravitational acceleration
vt=40*h0;           %Transition time in continuous friction function

%Application Data
    %app=1, Three Slider Mechanism
    %app=2, Slider-Crank
    %app=3, Quick Return Mechanism
    %app=4, Lumped Mass Coil Spring-5 masses
    %app=5,Rotating Disk with Translating Body
    %app=6, Double Pendulum
    %app=7, Loader
    
app=5;

[nb,ngc,nh,nc,nv,nu,NTSDA,NRSDA,PJDT,PMDT,PTSDAT,PRSDAT,q0,qd0]=...
    AppData(app);     %AppData(app)

par=[nb;ngc;nh;nc;nv;nu;g;utol;Btol;intol;Atol;Vtol;hvar;NTSDA;NRSDA;vt];

%Check initial conditions for feasibility
Phi0=PhiEval(0,q0,PJDT,par);
if norm(Phi0)>0.01
infeasableconfig=1;
return
end
Phiq0=PhiqEval(0,q0,PJDT,par);
Phid0=Phiq0*qd0;
if norm(Phid0)>0.01
infeasablevelocity=1; 
return
end

%Initial condition calculation, if required

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
Vvdd0=zeros(nv,10);
LLam0=zeros(nc,10);
Freac=zeros(nh,10);     %Vector of constraint reaction forces

Q(:,1)=q0;
Qd(:,1)=qd0;

%Evaluate Initial Acceleration and Lagrange Multipliers;

%Initial Parameterization
Phiq=PhiqEval(0,q0,PJDT,par);
U=Phiq';
B=inv(U'*U);
V=null(U');
v=zeros(nv,1);
vd=V'*qd0;
u=zeros(nu,1);
Vv(:,1)=v;
Vvd(:,1)=vd;
Uu(:,1)=u;

%Increment friction coefficients to obtain initial conditions on Lam and
%vdd
N=10;
w=0;
vdd=zeros(nv,1);
Lam=zeros(nc,1);
while w<=N
[vdd,Lam,jodeiter,ECond]=FrODEfunct0w(q0,qd0,vdd,Lam,...
    V,U,B,PJDT,PMDT,PTSDAT,PRSDAT,par,w,N);
w=w+1;
Econd0(w)=ECond;
jodeiter0(w)=jodeiter;
Vvdd0(:,w)=vdd;
LLam0(:,w)=Lam;

end

Gam=GamEval(0,q0,qd0,PJDT,par);
qdd=V*vdd-U*B*Gam;
Qdd(:,1)=qdd;
Vvdd(:,1)=vdd;
LLam(:,1)=Lam;

%Initialize Data For Integration
n=1;
t(1)=0;
vnorm(1)=Maxv+1;
jRepar=0;       %Counter for Reparameterization time step
Cr=2;           %Reparameterize if Cr>1
jiter=0;
Iterurpt(1)=0;
JCondrpt(1)=0;
ImpSoliterrpt(1)=MaxImpSoliter+1;
R1Normrpt(1)=R1nmax+1;
vnormrpt(1)=Maxv+1;
unormrpt(1)=Maxu+1;
ECondrpt(1)=MaxECond+1;
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
irpvnorm=0;
irpunorm=0;
irpECond=0;
irpBCond=0;
irpBCondRat=0;
irpvnorm=0;
irpunorm=0;
irpImpSoliterrpt=0;
irpJCond=0;
irpR1Norm=0;
irpBCond=0;

%Integration
while t(n)<tfinal
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);
tnm=t(n-1);

if integ<3     %Start implicit reparameterization criteria 
if BCondRat(n-1)>10
irpBCondRat=irpBCondRat+1;
Cr=2; 
end
if vnormrpt(n-1)>Maxv
irpvnorm=irpvnorm+1;
Cr=2;
end
if unormrpt(n-1)>Maxu
irpunorm=irpunorm+1;
Cr=2;
end
if ImpSoliterrpt(n-1)>MaxImpSoliter
irpImpSoliterrpt=irpImpSoliterrpt+1;
Cr=Cr+2;
end
if JCondrpt(n-1)>MaxJcond
irpJCond=irpJCond+1;
Cr=Cr+2;
end
if R1Normrpt(n-1)>R1nmax
irpR1Norm=irpR1Norm+1;
Cr=Cr+2;
end  
if BCond(n-1)>MaxBCond
irpBCond=irpBCond+1;
Cr=Cr+2;
end
end      %end implicit reparameterization criteria
if integ>2      %Start explicit reparameterization criteria    
if vnormrpt(n-1)>Maxv
irpvnorm=irpvnorm+1;
Cr=2;
end
if unormrpt(n-1)>Maxu
irpunorm=irpunorm+1;
Cr=2;
end    
if ECondrpt(n-1)>MaxECond
irpECond=irpECond+1;
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

if integ==1        %Implicit Index 0 trapezoidal
[v,vd,vdd,Lam,R1n,ImpSoliter,JCond,h,nch]=...
    ImplicitInd0Trap(n,tn,Vv,Vvd,Vvdd,Uu,LLam,q0,V,U,B,...
    h,hmax,nch,PMDT,PTSDAT,PRSDAT,PJDT,par);
    
R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
ImpSoliterrpt(n)=ImpSoliter;
hrpt(n)=h;    
end

if integ==2     %Implicit Index 0 SDIRK54
[v,vd,vdd,Lam,ImpSoliter,R1n,J,JCond,h,nch]=...
    ImplicitInd0SDIRK54(n,tn,Vv,Vvd,Vvdd,LLam,Uu,q0,U,V,B,par,...
    h,hmax,nch,npar,PJDT,PMDT,PTSDAT,PRSDAT);
R1Normrpt(n)=R1n;
JCondrpt(n)=JCond;
ImpSoliterrpt(n)=ImpSoliter;
hrpt(n)=h;
end

if integ==3         %Explicit Index0 Nystrom4
[v,vd,vdd,Lam,ECond,jodeiter1]=ExplicitInd0Nystrom4(n,tn,Vv,Vvd,...
    Vvdd,LLam,Uu,V,U,B,q0,h,npar,PJDT,PMDT,PTSDAT,PRSDAT,par); 
ECondrpt(n)=ECond;
jodeiter1rpt(n)=jodeiter1;
end

if integ==4         %Explicit Index0 RKFN45
[v,vd,vdd,Lam,ECond,Err,h,nch,jodeiter1]=...
    ExplicitInd0RKFN45(n,tn,Vv,Vvd,Vvdd,LLam,Uu,...
    U,V,B,q0,h,hmax,par,PMDT,PTSDAT,PRSDAT,PJDT,nch);
ECondrpt(n)=ECond;
Errrpt(n)=Err;
hrpt(n)=h;
jodeiter1rpt(n)=jodeiter1;
end

if integ==5         %Explicit Index0 Forward Euler
[v,vd,vdd,Lam,ECond,jodeiter]=ExplicitInd0ForwEuler(n,tn,Vv,...
    Vvd,Vvdd,LLam,Uu,V,U,B,q0,h,npar,PJDT,PMDT,PTSDAT,PRSDAT,par);

ECondrpt(n)=ECond;
jodeiterrpt(n)=jodeiter;
end

%Process Results
Vv(:,n)=v;
Vvd(:,n)=vd;
Vvdd(:,n)=vdd;
LLam(:,n)=Lam;
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

%Evaluate qdd
Gam=GamEval(tn,q,qd,PJDT,par);
qdd=D*vdd-U*B*Gam;
Qdd(:,n)=qdd;

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

if app==1   %Three Slider Mechanism
y1(n)=q(2); 
x2(n)=q(4); 
x3(n)=q(7);
y1d(n)=qd(2); 
x2d(n)=qd(4); 
x3d(n)=qd(7);
y1dd(n)=qdd(2); 
x2dd(n)=qdd(4); 
x3dd(n)=qdd(7);
x3m2(n)=q(7)-q(4);
P=[0,-1;1,0];
[i,j,s1pr,s2pr,v1pr,v2pr,mus,mud,d,ms,nm]=TranPart(1,PJDT);
[r1,ph1]=qPart(q,1);
A1=ATran(ph1);
F11(n)=Lam(1);
T11(n)=-s1pr'*v1pr*Lam(1)+[-r1'*A1*v1pr,v1pr'*A1'*v2pr]*[Lam(1);Lam(2)];
[i,j,s1pr,s2pr,v1pr,v2pr,mus,mud,d]=TranPart(2,PJDT);
[r1,ph1]=qPart(q,2);
A1=ATran(ph1);
F22(n)=Lam(3);
T22(n)=-s1pr'*v1pr*Lam(3)+[-r1'*A1*v1pr,v1pr'*A1'*v2pr]*[Lam(3);Lam(4)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Slider-Crank
phi1(n)=q(3);
phi1d(n)=qd(3);
phi1dd(n)=qdd(3);
x2(n)=q(4);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Quick Return
phi1(n)=q(3);
phi2(n)=q(6);
phi2d(n)=qd(6);
x4(n)=q(10);
x4d(n)=qd(10);
x4dd(n)=qdd(10);
[QA,F]=QAEval(tn,q,qd,Lam,PMDT,PJDT,PTSDAT,PRSDAT,par);
Freac(:,n)=F;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Lumped Mass Coil Spring 
delx1(n)=q(1)-0.2; 
delx2(n)=q(4)-0.4;
delx3(n)=q(7)-0.6;
delx4(n)=q(10)-0.8;
delx5(n)=q(13)-1;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5  %Rotating Disk with Translating Body
phi1rpt(n)=q(3);
A1=ATran(q(3));
r2=[q(4);q(5)];
dely2pr(n)=(A1'*(r2-A1*ux))'*uy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==6   %Double Pendulum
x1(n)=q(1); 
y1(n)=q(2); 
phi1(n)=q(3); 
phi1d(n)=qd(3);
x2(n)=q(4); 
y2(n)=q(5); 
phi2(n)=q(6);
phi2d(n)=qd(6);
phi2m1(n)=phi2(n)-phi1(n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==7   %Loader
y1(n)=q(2);
y2(n)=q(5);
phi1(n)=q(3);
phi2(n)=q(6);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate constraint error
Phi=PhiEval(tn,q,PJDT,par);
Phiq=PhiqEval(tn,q,PJDT,par);
Gam=GamEval(tn,q,qd,PJDT,par);
PosConstrNorm(n)=norm(Phi);
VelConstrNorm(n)=norm(Phiq*qd);
AccConstrNorm(n)=norm(Phiq*qdd+Gam);


end

