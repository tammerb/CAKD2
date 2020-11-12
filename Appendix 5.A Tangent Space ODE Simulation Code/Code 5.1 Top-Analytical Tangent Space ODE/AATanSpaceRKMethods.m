
%Code 5.1 Spinning Top With Analytical Tangent Space ODE

intol=10^-6;    %Tolerance in solving discretized equations of motion
Atol=10^-6;
vmax=0.75;      %Limit on norm of v%Absolute error tolerance for variable step metnods

h=0.001;        %Step size
hmax=0.01;     %Maximum Allowable Step size
hvar=2;         %hvar=1, variable step;hvar=2, fixed step

tfinal=100;     %Final time

Integ=1;    %Numerical Integration Options

    %Explicit Integration Methods:
        %Integ=1, Nystrom4; Integ=2, RungeKutta4;
        %Integ=3, Kutta3/8; integ=4, RKFN
        
nv=3;    % Variable dimension

% Problem Data
m=30;
g=9.8;
Jpr=diag([120;120;30]);
omegaz0=13.5;    %Initial angular velocity about z-axis
eps=10^-12;     %Angular velocity perturbation about x- and y-axes
Dcf=0;          %Drag coefficient
%Fixed Parameter Data List-Complete Partitioning in function parPart if
%explicit use of variable names is desired in user supplied functions. 
%Likewise variable partitions in qpart and qdpart

par=[nv;intol;Atol;hmax;hvar];

% Data Storage Arrays; Array named Vv, since V used in tangent space eqs.
Vv=zeros(nv,10);
Vvd=zeros(nv,10);
Vvdd=zeros(nv,10);
P=zeros(nv+1,10);
Pd=zeros(nv+1,10);
Pdd=zeros(nv+1,10);

% Initial Conditions
v0=[0;0;0];
p0=[1;0;0;0];
G0=Geval(p0);
pd0=0.5*G0'*[eps;eps;omegaz0];
V=G0';
vd0=V'*pd0;
Vv(:,1)=v0;
Vvd(:,1)=vd0;
P(:,1)=p0;
Pd(:,1)=pd0;

% Logic to reparameterize was REQUIRED BEYOND THIS POINT

%Integration Preparation
n=1;
t(1)=0;

%Initial Parameters
nch=1;      %Time of last change in step size
nrepar=0;   %Reparameterization counter

while t(n)<tfinal
%Time Step Update
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);

% Integration
tnm=t(n-1);
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);

%Test for reparameterization
if norm(vnm)>vmax
    p0=P(:,n-1);
    V=Geval(p0)';
    vnm=zeros(nv,1);
    vdnm=V'*Pd(:,n-1);
    Vv(:,n-1)=vnm;
    Vvd(:,n-1)=vdnm;
    nrepar=nrepar+1;
end
    
if Integ<5      %Explicit integrators
if Integ==1
[vn,vdn,vddn,Mcond]=ExplicitNystrom4(tnm,vnm,vdnm,h,par,p0,V,Jpr,m,g,Dcf);

end

if Integ==2
    [vn,vdn,vddn,Mcond]=ExplicitRungeKutta4(tnm,vnm,vdnm,h,par,...
        p0,V,Jpr,m,g,Dcf);
end

if Integ==3
    [vn,vdn,vddn,Mcond]=ExplicitKutta38(tnm,vnm,vdnm,h,par,...
        p0,V,Jpr,m,g,Dcf);
end

if Integ==4
[vn,vdn,vddn,Mcond,h,nch]=ExplicitRKFN45(n,tnm,vnm,vdnm,h,nch,...
    par,hvar,nv,Atol,hmax,p0,V,Jpr,m,g,Dcf);
hrpt(n)=h;
end

Mcondrpt(n)=Mcond;
end

%Evaluate and Record Solution
Vv(:,n)=vn;
vnorm(n)=norm(vn);
Vvd(:,n)=vdn;
vdnorm(n)=norm(vdn);
Vvdd(:,n)=vddn;
vddnorm(n)=norm(vddn);
pn=p0+V*vn-(1-sqrt(1-vn'*vn))*p0;
D=(eye(4)-(1/(pn'*p0))*p0*pn')*V;
pdn=D*vdn;
pddn=D*vddn-(pdn'*pdn/(pn'*p0))*p0;
P(:,n)=pn;
Pd(:,n)=pdn;
Pdd(:,n)=pddn;

%Report key data
omegpr=2*Geval(pn)*pdn;
omegaprnorm(n)=norm(omegpr);
uz=[0;0;1];
rn=BAeval(pn)*uz;
rx(n)=rn(1);
ry(n)=rn(2);
rz(n)=rn(3);
TE(n)=0.5*vdn'*AM1(vn,par,p0,V,Jpr,m,g)*vdn+m*g*rz(n);


end




   