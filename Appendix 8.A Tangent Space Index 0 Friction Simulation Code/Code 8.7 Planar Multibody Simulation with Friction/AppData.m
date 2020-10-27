function [nb,ngc,nh,nc,nv,nu,NTSDA,NRSDA,PJDT,PMDT,PTSDAT,...
    PRSDAT,q0,qd0]=AppData(app)

if app==1  %Three Body Translating Model
nb=3;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=2;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[2;1;0;z2;z2;0;0.1*ux;uy;0.1;0.5;0.3;1;2];    %Tran-Bod1 to ground
PJDT(:,2)=[2;2;0;z2;z2;0;0.1*ux;ux;0.1;0.5;0.3;3;2];   %Tran-Bod2 to ground
PJDT(:,3)=[2;3;0;z2;z2;0;0.1*ux;ux;0.1;0.5;0.3;5;2];   %Tran-Bod3 to ground
PJDT(:,4)=[3;1;2;z2;z2;5;z2;z2;zeros(3,1);7;1];       %Dist.-Bod1 to Bod2
nh=4;       %Number of holonomic constraints
nc=7;       %Number of constraint equations
nv=ngc-nc;
nu=nc;

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[6;1],[2;1],[2;1]];

%PTSDAT(10,NTSDA) TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT(:,1)=[1;0;z2;-10*uy;10;0;10;0];  %Bod1 to Grnd
PTSDAT(:,2)=[2;3;z2;10*ux;10;0;11;0];  %Bod2 to bod3

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0;0;pi/2;5*ux;0;6*ux;0];
qd0=[0;0;0;0;0;0;-1;0;0];

end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if app==2   %Slider-Crank

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=0;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[1;1;0;z2;z2;0;z2;z2;0.1;0.25;0.2;1;2];    %Revolute-crank to ground
PJDT(:,2)=[2;2;0;z2;z2;0;0.1*ux;ux;0.1;0.25;0.2;3;2];      %Trans.-slider2 to ground
PJDT(:,3)=[3;1;2;ux;z2;1.25;z2;z2;zeros(3,1);5;1];     %Dist.-crank to slider2
nh=3;       %Number of holonomic constraints
nc=5;      %Number of holonomic constraint equations
nv=ngc-nc;  %Number of independent coordinates
nu=nc;      %Number of dependent coordinates

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[5;5],[1;1]];

%PTSDAT(10,NTSDA) TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT=zeros(10,1);

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0;0;0;2.25*ux;0];
qd0=[0;0;100;0;0;0];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if app==3   %Quick Return

nb=4;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=0;    %Number of TSDA force elements
NRSDA=1;    %Number of RSDA torque elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[1;1;0;-2*ux;z2;0;z2;z2;0.1;0.3;0.25;1;2];    %Revolute-bar to ground
PJDT(:,2)=[1;2;0;z2;2*uy;0;z2;z2;0.1;0.3;0.25;3;2];    %Revolute-crank to ground
PJDT(:,3)=[1;2;3;1.5*ux;z2;0;z2;z2;0.1;0.3;0.25;5;2];    %Revolute-crank to key
PJDT(:,4)=[2;3;1;z2;z2;0;0.2*ux;ux;0;0.3;0.25;7;2];      %Trans.-bar to key
PJDT(:,5)=[2;4;0;z2;4*uy;0;0.1*ux;ux;0;0.3;0.25;9;2];    %Trans.-cutter to ground
PJDT(:,6)=[3;1;4;2*ux;z2;2.5298;z2;z2;zeros(3,1);11;1];  %Dist.-bar to cutter

nh=6;       %Number of holonomic constraints
nc=11;      %Number of holonomic constraint equations
nv=ngc-nc;  %Number of independent coordinates
nu=nc;      %Number of dependent coordinates

%PMDT(2,nb): Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[100;100],[1000;1000],[1;1],[50;50]];

%PTSDAT(10,NTSDA): TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies connected, si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=constant force 
PTSDAT=zeros(10,1);

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=RSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT(:,1)=[2;0;0;0;0;5000];

%Initial generalized coordinates
q10=[1.2;1.6;0.9273];
q20=[0;2;0];
q30=[1.5;2;0.9273];
q40=[0;4;0];

q0=[q10;q20;q30;q40];
qd0=zeros(12,1);  %Placeholder, qd0 calculated in main program

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Lumped Mass Coil Spring-5 masses
    
nb=5;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=5;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[2;1;0;z2;z2;0;0.1*ux;ux;0;0.2;0.19;1;2];      %Tran.-Body 1 to Ground
PJDT(:,2)=[2;2;0;z2;z2;0;0.1*ux;ux;0;0.2;0.19;3;2];      %Tran.-Body 2 to Ground
PJDT(:,3)=[2;3;0;z2;z2;0;0.1*ux;ux;0;0.2;0.19;5;2];      %Tran.-Body 3 to Ground
PJDT(:,4)=[2;4;0;z2;z2;0;0.1*ux;ux;0;0.2;0.19;7;2];      %Tran.-Body 4 to Ground
PJDT(:,5)=[2;5;0;z2;z2;0;0.1*ux;ux;0;0.2;0.19;9;2];      %Tran.-Body 5 to Ground

nh=5;       %Number of holonomic constraints
nc=10;       %Number of constraint equations
nv=ngc-nc;
nu=nc;

%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]] Mass Data TAble
PMDT(:,1)=[0.2;0.1];
PMDT(:,2)=[0.2;0.1];
PMDT(:,3)=[0.2;0.1];
PMDT(:,4)=[0.2;0.1];
PMDT(:,5)=[0.2;0.1];

%PTSDAT(10,NTSDA) Planar TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT(:,1)=[1;0;z2;z2;5000;0;0.2;0];    %Body 1 to Ground
PTSDAT(:,2)=[1;2;z2;z2;5000;0;0.2;0];    %Body 1 to 2
PTSDAT(:,3)=[2;3;z2;z2;5000;0;0.2;0];    %Body 2 to 3
PTSDAT(:,4)=[3;4;z2;z2;5000;0;0.2;0];    %Body 3 to 4
PTSDAT(:,5)=[4;5;z2;z2;5000;0;0.2;0];    %Body 4 to 5

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0.2;0;0;0.4;0;0;0.6;0;0;0.8;0;0;1;0;0];
qd0=[zeros(12,1);-1;0;0];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if app==5   %Rotating Disk with Translating Body

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=1;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[1;1;0;z2;z2;0;z2;z2;0.1;0.05;0.04;1;2];    %Revolute-1 to ground
PJDT(:,2)=[2;1;2;ux;z2;0;0.1*uy;uy;0;0.05;0.04;3;2];   %Trans.-1 to 2

nh=2;       %Number of holonomic constraints
nc=4;      %Number of holonomic constraint equations
nv=ngc-nc;  %Number of independent coordinates
nu=nc;      %Number of dependent coordinates

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[10;10],[5;5]];

%PTSDAT(10,NTSDA) TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT(:,1)=[1;2;ux+10*uy;z2;10;0;10.1;0];

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0;0;0;ux;0];
qd0=[0;0;0;0;0;0];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if app==6       %Double Pendulum

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=1;    %Number of TSDA force elements
NRSDA=2;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[1;1;0;-ux;z2;0;z2;z2;0.1;0.2;0.1;1;2];     %Rev-Body 1 to Ground
PJDT(:,2)=[1;1;2;ux;-ux;0;z2;z2;0.1;0.2;0.1;3;2];     %Rev-Body 1 to Body 2

nh=2;       %Number of holonomic constraints
nc=4;       %Number of constraint equations
nv=ngc-nc;
nu=nc;

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[1;0.3],[1;0.3]];

%PTSDAT(10,NTSDA) TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
%Unilateral spring is defined by setting spring constant in QA, QAsq,
%QAsdd and TE calculation equal to K*(1-sign(q(1)))/2
K=0;
PTSDAT(:,1)=[1;0;z2;-ux-uy;K;0;1;0]; 

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT(:,1)=[1,0,0,0,pi/2,0];
PRSDAT(:,2)=[1,2,0,0,0,0];

%Initial generalized coordinates

q0=[1;0;0;3;0;0];
qd0=[0;1;1;0;3;1];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if app==7   %Loader
    
nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
NTSDA=4;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
z2=zeros(2,1);

%PJDT(17,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr;R;mus;mud;ms;nm]; 
    %k=joint No., t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis,
    %length of vi is di in Eq. (6.4.17), R=rad Rev, mus&mud=FrCoefs, 
    %ms=Lagrange muplt.start address, nm=no. of mujlt.
PJDT(:,1)=[1;1;2;uy;3*ux+1.5*uy;0;z2;z2;zeros(3,1);1;2];      %Rev.-Body to Loader
nh=1;       %Number of holonomic constraints
nc=2;      %Number of holonomic constraint equations
nv=ngc-nc;  %Number of independent coordinates
nu=nc;      %Number of dependent coordinates

%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]] Mass Data TAble
PMDT(:,1)=[2000;2000];
PMDT(:,2)=[750;100];
 
%PTSDAT(10,NTSDA) Planar TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT(:,1)=[1;0;ux-uy;ux-5*uy;3*10^5;1000;4;0];
PTSDAT(:,2)=[1;0;ux-uy;5*ux-uy;3*10^5;1000;4;0];
PTSDAT(:,3)=[1;0;-2*ux-uy;-2*ux-5*uy;6*10^5;1000;4;0];
PTSDAT(:,4)=[1;2;0.5*uy;2*ux+uy;0;0;1;-50000];

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0;0;0;-3*ux-0.5*uy;0];
qd0=[0;0;0;0;0;0];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end









