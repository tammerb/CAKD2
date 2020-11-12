function [nb,ngc,nh,nc,NTSDA,NRSDA,PJDT,PMDT,PTSDAT,PRSDAT,q0,qd0]=...
    AppData(app)

if app==1   %Quick Return

nb=4;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=6;       %Number of holonomic constraints
nc=11;      %Number of constraint equations
nv=ngc-nc;  %Number of independent coordinates
nu=nc;      %Number of dependent coordinates
NTSDA=0;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh): Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;-2*ux;zer;0;zer;zer];    %Revolute-bar to ground
PJDT(:,2)=[1;2;0;zer;2*uy;0;zer;zer];    %Revolute-crank to ground
PJDT(:,3)=[1;2;3;1.5*ux;zer;0;zer;zer];    %Revolute-crank to key
PJDT(:,4)=[2;1;3;zer;zer;0;ux;ux];      %Trans.-bar to key
PJDT(:,5)=[2;4;0;zer;4*uy;0;ux;ux];      %Trans.-cutter to ground
PJDT(:,6)=[3;1;4;2*ux;zer;2.5298;zer;zer];     %Dist.-bar to cutter

%PMDT(2,nb): Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[100;100],[1000;1000],[1;1],[50;50]];

%PTSDAT(10,NTSDA): TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies connected, si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=constant force 
PTSDAT=zeros(10,1);

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates
q10=[1.2;1.6;0.9273];
q20=[0;2;0];
q30=[1.5;2;0.9273];
q40=[0;4;0];

q0=[q10;q20;q30;q40];
qd0=zeros(12,1);  %Placeholder, qd0 calculated in main program

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Loader
    
nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=1;       %Number of holonomic constraints
nc=2;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=4;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist,4=RevClr), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;2;uy;3*ux+1.5*uy;0;zer;zer];      %Rev.-Body to Loader

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Slider-Crank

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nc=5;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;zer;zer;0;zer;zer];    %Revolute-crank to ground
PJDT(:,2)=[2;2;0;zer;zer;0;ux;ux];      %Trans.-slider2 to ground
PJDT(:,3)=[3;1;2;ux;zer;2;zer;zer];     %Dist.-crank to slider2

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

q0=[0;0;0;3*ux;0];
qd0=[0;0;100;0;0;0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %multiple slider-crank
    
nb=6;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=11;       %Number of holonomic constraints
nc=17;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis

PJDT(:,1)=[1;1;0;zer;zer;0;zer;zer];    %Revolute-crank to ground
PJDT(:,2)=[2;2;0;zer;zer;0;ux;ux];      %Trans.-slider2 to ground
PJDT(:,3)=[2;3;0;zer;zer;0;ux;ux];      %Trans.-slider3 to ground
PJDT(:,4)=[2;4;0;zer;zer;0;ux;ux];      %Trans.-slider4 to ground
PJDT(:,5)=[2;5;0;zer;zer;0;ux;ux];      %Trans.-slider5 to ground
PJDT(:,6)=[2;6;0;zer;zer;0;ux;ux];      %Trans.-slider6 to ground
PJDT(:,7)=[3;1;2;ux;zer;2;zer;zer];     %Dist.-crank to slider2
PJDT(:,8)=[3;1;3;ux;zer;2;zer;zer];     %Dist.-crank to slider3
PJDT(:,9)=[3;1;4;ux;zer;2;zer;zer];     %Dist.-crank to slider4
PJDT(:,10)=[3;1;5;ux;zer;2;zer;zer];     %Dist.-crank to slider5
PJDT(:,11)=[3;1;6;ux;zer;2;zer;zer];     %Dist.-crank to slider6

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[5;5],[1;1],[1;1],[1;1],[1;1],[1;1]];

%PTSDAT(10,NTSDA) TSDA Data Tabl
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

q0=zeros(ngc,1);
i=2;
while i<=6
    qs=[3*ux;0];
    q0=Add(q0,qs,3*(i-1),0);
    i=i+1;
end
qd0=[0;0;100;zeros(15,1)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5       %Pendulum-Distance constraint

nb=1;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=1;       %Number of holonomic constraints
nc=1;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[3;1;0;zer;zer;1;zer;zer];     %Dist.-Body CG to Ground

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[1;1]];

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

q0=[1;0;0];
qd0=[0;0;0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==6       %Pendulum-Revolute constraint

nb=1;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=1;       %Number of holonomic constraints
nc=2;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;-ux;zer;0;zer;zer];     %Rev-Body 1 to Ground

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[1;0.001]];

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

q0=[1;0;0];
qd0=[0;0;0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==7       %Double Pendulum-Unilateral Spring

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=2;       %Number of holonomic constraints
nc=4;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=1;    %Number of TSDA force elements
NRSDA=2;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;-ux;zer;0;zer;zer];     %Rev-Body 1 to Ground
PJDT(:,2)=[1;1;2;ux;-ux;0;zer;zer];     %Rev-Body 1 to Body 2

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[1;0.3],[1;0.3]];

%PTSDAT(10,NTSDA) TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
%Unilateral spring is defined by setting spring constant in QA, QAsq,
%and QAsdd equal to K=K*(1-sign(q(1)))/2
K3=10^5;
PTSDAT(:,1)=[1;0;zer;-ux-uy;K3;0;1;0]; 

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT(:,1)=[1,0,0,0,pi/2,0];
PRSDAT(:,2)=[1,2,10^5,10^4,0,0];

%Initial generalized coordinates

q0=[1;0;0;3;0;0];
qd0=[0;0;0;0;0;0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==8       %Double Pendulum

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=2;       %Number of holonomic constraints
nc=4;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements
NRSDA=2;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;-ux;zer;0;zer;zer];     %Rev-Body 1 to Ground
PJDT(:,2)=[1;1;2;ux;-ux;0;zer;zer];     %Rev-Body 1 to Body 2

%PMDT(2,nb) Mass Data Table
%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]]
PMDT=[[1;0.3],[1;0.3]];

%PTSDAT(10,NTSDA) TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT=zeros(10,1); 

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT(:,1)=[1,0,20,0,pi/2,0];
PRSDAT(:,2)=[1,2,10^3,0,0,0];

%Initial generalized coordinates

q0=[0;-1;-pi/2;0;-3;-pi/2];
qd0=[0;0;0;10;0;10];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==9   %Lumped Mass Coil Spring-5 masses
    
nb=5;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=5;       %Number of holonomic constraints
nc=10;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=5;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist,4=RevClr), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[2;1;0;zer;zer;0;ux;ux];      %Tran.-Body 1 to Ground
PJDT(:,2)=[2;1;2;zer;zer;0;ux;ux];      %Tran.-Body 1 to 2
PJDT(:,3)=[2;2;3;zer;zer;0;ux;ux];      %Tran.-Body 2 to 3
PJDT(:,4)=[2;3;4;zer;zer;0;ux;ux];      %Tran.-Body 3 to 4
PJDT(:,5)=[2;4;5;zer;zer;0;ux;ux];      %Tran.-Body 4 to 5

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
PTSDAT(:,1)=[1;0;zer;zer;5000;0;0.2;0];    %Body 1 to Ground
PTSDAT(:,2)=[1;2;zer;zer;5000;0;0.2;0];    %Body 1 to 2
PTSDAT(:,3)=[2;3;zer;zer;5000;0;0.2;0];    %Body 2 to 3
PTSDAT(:,4)=[3;4;zer;zer;5000;0;0.2;0];    %Body 3 to 4
PTSDAT(:,5)=[4;5;zer;zer;5000;0;0.2;0];    %Body 4 to 5

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0.2;0;0;0.4;0;0;0.6;0;0;0.8;0;0;1;0;0];
qd0=[zeros(12,1);-1;0;0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==10   %Lumped Mass Coil Spring-10 masses
    
nb=10;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=10;       %Number of holonomic constraints
nc=20;       %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=10;    %Number of TSDA force elements
NRSDA=0;    %Number of RSDA force elements

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh) Joint Data Table
%PJTd(:,k)=[t;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %t=joint type(1=Rev,2=Tran,3=Dist,4=RevClr), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[2;1;0;zer;zer;0;ux;ux];      %Tran.-Body 1 to Ground
PJDT(:,2)=[2;1;2;zer;zer;0;ux;ux];      %Tran.-Body 1 to 2
PJDT(:,3)=[2;2;3;zer;zer;0;ux;ux];      %Tran.-Body 2 to 3
PJDT(:,4)=[2;3;4;zer;zer;0;ux;ux];      %Tran.-Body 3 to 4
PJDT(:,5)=[2;4;5;zer;zer;0;ux;ux];      %Tran.-Body 4 to 5
PJDT(:,6)=[2;5;6;zer;zer;0;ux;ux];      %Tran.-Body 5 to 6
PJDT(:,7)=[2;6;7;zer;zer;0;ux;ux];      %Tran.-Body 6 to 7
PJDT(:,8)=[2;7;8;zer;zer;0;ux;ux];      %Tran.-Body 7 to 8
PJDT(:,9)=[2;8;9;zer;zer;0;ux;ux];      %Tran.-Body 8 to 9
PJDT(:,10)=[2;9;10;zer;zer;0;ux;ux];      %Tran.-Body 9 to 10


%PMDT=[[m1;J1],[m2;J2],...,[mnb;Jnb]] Mass Data TAble
PMDT(:,1)=[0.1;0.1];
PMDT(:,2)=[0.1;0.1];
PMDT(:,3)=[0.1;0.1];
PMDT(:,4)=[0.1;0.1];
PMDT(:,5)=[0.1;0.1];
PMDT(:,6)=[0.1;0.1];
PMDT(:,7)=[0.1;0.1];
PMDT(:,8)=[0.1;0.1];
PMDT(:,9)=[0.1;0.1];
PMDT(:,10)=[0.1;0.1];

%PTSDAT(10,NTSDA) Planar TSDA Data Table
%PTSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F]; T=TSDA No., 
    %i&j=bodies conn.,si&jpr=vectors to Pi&j, K=spring constant,
    %C=damping coefficient,el0=spring free length,F=const. force
PTSDAT(:,1)=[1;0;zer;zer;10000;0;0.1;0];    %Body 1 to Ground
PTSDAT(:,2)=[1;2;zer;zer;10000;0;0.1;0];    %Body 1 to 2
PTSDAT(:,3)=[2;3;zer;zer;10000;0;0.1;0];    %Body 2 to 3
PTSDAT(:,4)=[3;4;zer;zer;10000;0;0.1;0];    %Body 3 to 4
PTSDAT(:,5)=[4;5;zer;zer;10000;0;0.1;0];    %Body 4 to 5
PTSDAT(:,6)=[5;6;zer;zer;10000;0;0.1;0];    %Body 5 to 6
PTSDAT(:,7)=[6;7;zer;zer;10000;0;0.1;0];    %Body 6 to 7
PTSDAT(:,8)=[7;8;zer;zer;10000;0;0.1;0];    %Body 7 to 8
PTSDAT(:,9)=[8;9;zer;zer;10000;0;0.1;0];    %Body 8 to 9
PTSDAT(:,10)=[9;10;zer;zer;10000;0;0.1;0];    %Body 9 to 10

%PRSDAT(6,NRSDA): RSDA Data Table
%PRSDAT(:,R)=[i;j;K;C;phi0;T]; R=TSDA No., 
    %i&j=bodies connected, K=spring constant,
    %C=damping coefficient,phi0=spring free angle,T=constant torque 
PRSDAT=zeros(6,1);

%Initial generalized coordinates

q0=[0.1;0;0;0.2;0;0;0.3;0;0;0.4;0;0;0.5;0;0;0.6;0;0;0.7;0;0;0.8;0;0;...
    0.9;0;0;1;0;0];
qd0=[zeros(27,1);-1;0;0];

end

end







