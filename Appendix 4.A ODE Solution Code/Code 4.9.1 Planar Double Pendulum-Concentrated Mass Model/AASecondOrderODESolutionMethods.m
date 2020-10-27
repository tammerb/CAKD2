%4.9.1 Planar Double Pendulum; Based on Code 4.8, Second Order ODE Solver

intol=10^-6;    %Tolerance in solving discretized equations of motion
Atol=10^-6;     %Absolute error tolerance for variable step metnods

h=0.0001;        %Step size
hmax=0.001;     %Maximum Allowable Step size
hvar=1;         %hvar=1, variable step;hvar=2, fixed step

tfinal=9.99;     %Final time

Integ=6;    %Numerical Integration Options

    %Explicit Integration Methods:
        %Integ=1, Nystrom4; Integ=2, RungeKutta4;
        %Integ=4, Kutta3/8; integ=4, RKFN
                                    
    %Implicit Integration Methods:
        %Integ=5, Trapezoidal; Integ=6, SDIRK54
        
nv=2;      %User define variable dimension

%Fixed Parameter Data List-Complete Partitioning in function parPart if
%explicit use of variable names is desired in user supplied functions. 
par=[nv;intol;Atol;hmax;hvar];

% Problem Data
m1=1;
m2=100;

%Problem Parameter Data List-Partitioning in function AdatPart if
%use of variable names is desired in user supplied functions.
dat=[m1;m2];

% Data Storage Arrays
V=zeros(nv,10);
Vd=zeros(nv,10);
Vdd=zeros(nv,10);

% Initial Conditions
v0=[0;0];   %User define initial position
vd0=[0;0];  %User define initial velocity

V(:,1)=v0;
Vd(:,1)=vd0;

% NO USER CHANGES/INPUT REQUIRED BEYOND THIS POINT

%Integration Preparation
n=1;
t(1)=0;

%Initial Parameters
nch=1;

% Integration

while t(n)<tfinal
%Time Step Update
n=n+1;
t(n)=t(n-1)+h;
tn=t(n);


% Integration
tnm=t(n-1);
vnm=V(:,n-1);
vdnm=Vd(:,n-1);


if Integ<5      %Explicit integrators
if Integ==1
[vn,vdn,vddn,Mcond]=ExplicitNystrom4(tnm,vnm,vdnm,h,par,dat);

end

if Integ==2
    [vn,vdn,vddn,Mcond]=ExplicitRungeKutta4(tnm,vnm,vdnm,h,par,dat);
end

if Integ==3
    [vn,vdn,vddn,Mcond]=ExplicitKutta38(tnm,vnm,vdnm,h,par,dat);
end

if Integ==4
[vn,vdn,vddn,Mcond,h,nch]=ExplicitRKFN45(n,tnm,vnm,vdnm,h,nch,...
    par,dat,hvar,nv,Atol,hmax);
hrpt(n)=h;
end

Mcondrpt(n)=Mcond;
end

if Integ >4     %Implicit integrators
    
if Integ==5
[vn,vdn,vddn,jiter,R1Norm,JCond]=ImplicitTrap(n,tn,...
    V,Vd,Vdd,intol,par,dat,h);

jiterrpt(n)=jiter;
        
end
if Integ==6
    
[vn,vdn,vddn,Maxjiter,R1Norm,JCond,h,err]=ImplicitSDIRK54(n,tn,...
    V,Vd,Vdd,par,dat,intol,Atol,nv,h,hmax,nch,hvar); 

hrpt(n)=h;
R1Normrpt(n)=R1Norm;
errrpt(n)=err;
jiterrpt(n)=Maxjiter;
    
end

JCondrpt(n)=JCond;


end

%Evaluate and Record Solution
V(:,n)=vn;
vnorm(n)=norm(vn);
Vd(:,n)=vdn;
vdnorm(n)=norm(vdn);
Vdd(:,n)=vddn;
vddnorm(n)=norm(vddn);



%Report key data
v1(n)=vn(1);
v2(n)=vn(2);
v1d(n)=vdn(1);
v2d(n)=vdn(2);
v1dd(n)=vddn(1);
v2dd(n)=vddn(2);
s1=sin(vn(1));
s2=sin(vn(2));

TE(n)=0.5*vdn'*AM(vn,par,dat)*vdn+9.8*((m1+m2)*s1+m2*s2);


end




   