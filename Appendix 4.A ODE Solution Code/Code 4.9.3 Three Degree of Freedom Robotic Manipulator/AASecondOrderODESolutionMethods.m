
%4.7.3 Three Degree of Freedom Robotic Manipulator-Matrix Formulation, Eq.
%(4.9.6)

intol=10^-6;    %Tolerance in solving discretized equations of motion
Atol=10^-6;     %Absolute error tolerance for variable step metnods

h=0.001;        %Step size
hmax=0.01;     %Maximum Allowable Step size
hvar=2;         %hvar=1, variable step;hvar=2, fixed step

tfinal=10;     %Final time

integ=6;    %Numerical Integration Options

    %Explicit Integration Methods:
        %Integ=1, Nystrom4; Integ=2, RungeKutta4;
        %Integ=3, Kutta3/8; integ=4, RKFN
                                    
    %Implicit Integration Methods:
        %Integ=5, Trapezoidal; Integ=6, SDIRK54
        
nv=3;    % Variable dimension

%Fixed Parameter Data List-Complete Partitioning in function parPart if
%explicit use of variable names is desired in user supplied functions. 
%Likewise variable partitions in qpart and qdpart

par=[nv;intol;Atol;hmax;hvar];

% Problem Data
g=9.8;
m1=200;
m2=100;
m3=100;
%J1 is J1ii*I3
J1ii=40;
%J3=J2, with diagonal elements J21, J22, and J23
J21=20;
J22=5;
J23=20;
K1=1000;
K2=10000;
K3=10000;

%Problem Parameter Data List-Partitioning in function AdatPart if
%use of variable names is desired in user supplied functions.
dat=[g;m1;m2;m3;J1ii;J21;J22;J23;K1;K2;K3];

% Data Storage Arrays
V=zeros(nv,10);
Vd=zeros(nv,10);
Vdd=zeros(nv,10);

% Initial Conditions
v0=[0;0;0];
vd0=[1;0;0];
V(:,1)=v0;
Vd(:,1)=vd0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


if integ<5      %Explicit integrators
if integ==1
[vn,vdn,vddn,Mcond]=ExplicitNystrom4(tnm,vnm,vdnm,h,par,dat,integ);

end

if integ==2
    [vn,vdn,vddn,Mcond]=ExplicitRungeKutta4(tnm,vnm,vdnm,h,par,dat,integ);
end

if integ==3
    [vn,vdn,vddn,Mcond]=ExplicitKutta38(tnm,vnm,vdnm,h,par,dat,integ);
end

if integ==4
[vn,vdn,vddn,Mcond,h,nch]=ExplicitRKFN45(n,tnm,vnm,vdnm,h,nch,...
    par,dat,hvar,nv,Atol,hmax,integ);
hrpt(n)=h;
end

Mcondrpt(n)=Mcond;
end

if integ >4     %Implicit integrators
    
if integ==5
[vn,vdn,vddn,jiter,R1Norm,JCond]=ImplicitTrap(n,tn,...
    V,Vd,Vdd,intol,par,dat,h,integ);

jiterrpt(n)=jiter;
        
end
if integ==6
    
[vn,vdn,vddn,Maxjiter,R1Norm,JCond,h,nch,k1]=ImplicitSDIRK54(n,...
    tn,V,Vd,Vdd,par,dat,intol,Atol,nv,h,hmax,hvar,integ,nch); 

hrpt(n)=h;
R1Normrpt(n)=R1Norm;
jiterrpt(n)=Maxjiter;
k1mvddNorm(n)=norm(k1-vddn);
    
end

JCondrpt(n)=JCond;


end

%Record Solution
V(:,n)=vn;
vnorm(n)=norm(vn);
Vd(:,n)=vdn;
vdnorm(n)=norm(vdn);
Vdd(:,n)=vddn;
vddnorm(n)=norm(vddn);



%Report key data
v1(n)=vn(1);
v2(n)=vn(2);
v3(n)=vn(3);
v1d(n)=vdn(1);
v2d(n)=vdn(2);
v3d(n)=vdn(3);
v1dd(n)=vddn(1);
v2dd(n)=vddn(2);
v3dd(n)=vddn(3);
s1=sin(vn(1));
s2=sin(vn(2));
c1=cos(vn(1));
c2=cos(vn(2));
deriv=1;
[M,gf,M2,gfsv,gfsvd] = AMg(tn,vn,vdn,vddn,par,dat,integ,deriv);
KE(n)=0.5*vdn'*M*vdn;
PE(n)=m2*g*(2+s2)+m3*g*(2+(v3(n)+2)*s2)...
    +0.5*(K1*v1(n)^2+K2*v2(n)^2+K3*v3(n)^2);
TE(n)=KE(n)+PE(n);


end




   