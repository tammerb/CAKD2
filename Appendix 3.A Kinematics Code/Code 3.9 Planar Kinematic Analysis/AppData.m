function [nb,ngc,nh,nhc,nd,PJDT,q0e]=AppData(app)

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==1   %Slider-Crank

nb=2;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=5;       %Number of constraint equations
nd=ngc-nhc;  %Number of drivers

%PJDT(12,nh): Planar Joint Data Table (First nh joints not time dependent)
%PJDT(:,k)=[T;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
%T=joint type(1=Rev,2=Tran,3=Dist, 4=RotD), i&j=bodies connected,
%sipr&sjpr=vectors to Pi&Pj in joint definition, d=dist., 
%vipr&vjpr=vectors along translstional axis
PJDT(:,1)=[1;1;0;zer;zer;0;zer;zer];%Revolute-crank to ground
PJDT(:,2)=[2;2;0;zer;zer;0;ux;ux];%Trans.-slider2 to ground
PJDT(:,3)=[3;1;2;ux;zer;1.01;zer;zer];%Dist.-crank to slider2
PJDT(:,4)=[4;1;0;zer;zer;0;zer;zer];%RelRotDrver-crank to ground

q0e=[0;0;0;2.2*ux;0];  %Initial generalized coordinate estimate

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Quick Return

nb=4;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=6;       %Number of time independent holonomic constraints
nhc=11;     %Number of time independent holonomic constraint equations
nd=ngc-nhc; %Number of time dependent driving constraint equations


ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh): Planar Joint Data Table (First nh not time dependent)
%PJTd(:,k)=[T;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %T=joint type(1=Rev,2=Tran,3=Dist, 4=RotD), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;-2*ux;zer;0;zer;zer];    %Revolute-bar to ground
PJDT(:,2)=[1;2;0;zer;2*uy;0;zer;zer];    %Revolute-crank to ground
PJDT(:,3)=[1;2;3;1.5*ux;zer;0;zer;zer];    %Revolute-crank to key
PJDT(:,4)=[2;1;3;zer;zer;0;ux;ux];      %Trans.-bar to key
PJDT(:,5)=[2;4;0;zer;4*uy;0;ux;ux];      %Trans.-cutter to ground
PJDT(:,6)=[3;1;4;2*ux;zer;2.5298;zer;zer];     %Dist.-bar to cutter
PJDT(:,7)=[4;2;0;zer;zer;0;zer;zer];    %RotD-crank to ground

%Initial generalized coordinate estimate
q10=[1.2;1.6;0.9273];
q20=[0;2;0];
q30=[1.5;2;0.9273];
q40=[0;4;0];
q0e=[q10;q20;q30;q40];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Windshield Wiper

nb=3;       %Number of bodies
ngc=3*nb;   %number of generalized coordinates
nh=5;       %Number of time independent holonomic constraints
nhc=8;     %Number of time independent holonomic constraint equations
nd=ngc-nhc; %Number of time dependent driving constraint equations

ux=[1;0];
uy=[0;1];
zer=zeros(2,1);

%PJDT(12,nh): Planar Joint Data Table (First nh not time dependent)
%PJTd(:,k)=[T;i;j;sipr;sjpr;d;vipr;vjpr]; k=joint No., 
    %T=joint type(1=Rev,2=Tran,3=Dist, 4=RotD), i&j=bodies conn.,
    %si&jpr=vectors to Pi&j, d=dist., vi&jpr=vectors along trans axis
PJDT(:,1)=[1;1;0;zer;zer;0;zer;zer];    %Rev-bod 1 to ground
PJDT(:,2)=[1;2;0;zer;-70*ux;0;zer;zer]; %Rev-bod 2 to ground
PJDT(:,3)=[1;3;0;zer;30*ux;0;zer;zer];  %Rev-bod 3 to ground
PJDT(:,4)=[3;1;2;-5*uy;-8*uy;70;zer;zer]; %Dist.-bod 1 to bod 2
PJDT(:,5)=[3;2;3;5*uy;5*uy;100;zer;zer]; %Dist.-bod 2 to bod 3
PJDT(:,6)=[4;1;0;zer;zer;0;zer;zer];    %RotD-bod 1 to ground

%Initial generalized coordinate estimate
q10=[0;0;0];
q20=[-70;0;0];
q30=[30;0;0];
q0e=[q10;q20;q30];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end







