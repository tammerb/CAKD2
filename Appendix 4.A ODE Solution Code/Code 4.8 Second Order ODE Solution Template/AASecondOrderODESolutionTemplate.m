%Code 4.8: Solution of Second Order ODE M(v)vdd=g(v,vd,t) of Eq.(4.8.43)
%on time interval [0,tf]

%Functions beginning with A require user input; Functions beginning with
%B, E, I, and O do not require user input

%Function g=Agf(t,v,vd,par) defines g(v,vd,t); required for all integrators

%Function M=AM(v,par) defines M(v); required for all integrators

%Function [gsv,gsvd] = Agfsvvd(t,v,vd,par) defines partial derivatives of
%g(v,vd,t) w.r.t. v and vd; required for implicit integrators

%Function M2=AM2(v,mu,par) defines derivative of (M(v)mu) w.r.t. v;
%required for implicit integrators

%Function [,,,,]=AdatPart(dat) defines problem dependent data; may be
%defined if user wishes to employ named problem data in subroutines 

intol=10^-6;    %Tolerance in solving discretized equations
Atol=10^-6;     %Absolute error tolerance for variable step metnods

h=0.001;        %Step size
hmax=0.01;     %Maximum Allowable Step size
hvar=2;         %hvar=1, variable step;hvar=2, fixed step

tfinal=12;     %Final time

Integ=6;    %Numerical Integration Options

    %Explicit Integration Methods:
        %Integ=1, Nystrom4; Integ=2, RungeKutta4;
        %Integ=3, Kutta3/8; integ=4, RKFN
                                    
    %Implicit Integration Methods:
        %Integ=5, Trapezoidal; Integ=6, SDIRK54
        
nv=[];    %User define variable dimension

%Fixed Parameter Data List-Partitioning in function BparPart if
%use of variable names is desired in user supplied functions. 
par=[nv;intol;Atol;hmax;hvar];

% Problem Data
%enter problem data by "name=value;" and enter into dat list.

%Problem Parameter Data List-Partitioning in function AdatPart if
%use of variable names is desired in user supplied functions. 
dat=[;;;;]; %User input data parameters here and in function AdatPart,
            %as in function BparPart
%Example:  dat=[m1,J1,k1]

% Data Storage Arrays
V=zeros(nv,10);
Vd=zeros(nv,10);
Vdd=zeros(nv,10);

% Initial Conditions
v0=[];      %User define initial position
vd0=[];     %User define initial velocity

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
s1=sin(vn(1));
s2=sin(vn(2));
c1=cos(vn(1));
c2=cos(vn(2));
TE(n)=0.5*vdn'*AM(vn,par,dat)*vdn+9.8*(2*s1+s2);


end




   