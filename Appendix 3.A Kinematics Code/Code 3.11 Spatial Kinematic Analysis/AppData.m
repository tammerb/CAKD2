function [nb,ngc,nh,nhc,nc,nd,SJDT,q0e]=AppData(app)

z3=[0;0;0];
z4=[0;0;0;0];
I3=eye(3);
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==1   %2 Bar with 2 rotational drivers
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=2;       %Number of holonomic constraints
nhc=10;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nd=ngc-nc;

%SJDT(22,nh): Spatial Joint Data Table 
%SJDT(:,k)=[T;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; T=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph, 9=DistDr, 10=RotDr); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[4;1;0;z3;z3;0;ux;uz;ux;uz];   %Rev-Body1 to Ground
SJDT(:,2)=[4;1;2;1.5*ux;z3;0;ux;uz;ux;uz];              %Tran-Body2 to Ground
SJDT(:,3)=[10;1;0;z3;z3;0;ux;uz;ux;uz];  %RotDr-Body 1 to ground
SJDT(:,4)=[10;1;2;1.5*ux;z3;0;ux;uz;ux;uz];  %RotDr-Body 1 to 2

%Initial generalized coordinate estimate
r10=z3;
p10=[1;0;0;0];
r20=1.5*ux;
p20=[1;0;0;0];
q0e=[r10;p10;r20;p20];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %4 Bar with internal rotational driver
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nd=ngc-nc;

%SJDT(22,nh): Spatial Joint Data Table 
%SJTd(:,k)=[T;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; T=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph, 9=DistDr, 10=RotDr); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[4;1;0;z3;z3;0;ux;uz;ux;uz];   %Rev-Body1 to Ground
SJDT(:,2)=[4;1;2;1.5*ux;z3;0;ux;uz;ux;uz];              %Tran-Body2 to Ground
SJDT(:,3)=[1;2;0;1.5*ux;4*uy;5;z3;z3;z3;z3];  %Dist-Body 1 to 2
SJDT(:,4)=[10;1;2;1.5*ux;z3;0;ux;uz;ux;uz];  %RotDr-Body 1 to 2

%Initial generalized coordinate estimate
r10=z3;
p10=[1;0;0;0];
r20=1.5*ux;
p20=[1;0;0;0];
q0e=[r10;p10;r20;p20];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %2-body Slider, x-z plane, x1 distDr
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nd=ngc-nc;

%SJDT(22,nh): Spatial Joint Data Table 
%SJTd(:,k)=[T;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; T=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph, 9=DistDr, 10=RotDr); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[5;1;0;z3;z3;0;uy;ux;uy;ux];   %Tran-Body1 to Ground
SJDT(:,2)=[5;2;0;z3;z3;0;ux;uz;ux;uz];    %Tran-Body2 to Ground
SJDT(:,3)=[1;1;2;z3;z3;5;z3;z3;z3;z3];  %Dist-Body 1 to 2
SJDT(:,4)=[9;1;0;z3;z3;0;z3;z3;z3;z3];  %DistDr-Body 1 to Ground

%Initial generalized coordinate estimate
r10=[3;0;0];
p10=[1;0;0;0];
r20=[0;0;4];
p20=[1;0;0;0];
q0e=[r10;p10;r20;p20];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Slider-Crank in y-z Plane, Bod 1 RotDr about x axis
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nd=ngc-nc;

%SJDT(22,nh): Spatial Joint Data Table 
%SJTd(:,k)=[T;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; T=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph, 9=DistDr, 10=RotDr); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[4;1;0;z3;z3;0;uy;ux;uy;ux];   %Rev-Body1 to Ground
SJDT(:,2)=[5;2;0;z3;z3;0;uz;uy;uz;uy];              %Tran-Body2 to Ground
SJDT(:,3)=[1;1;2;uy;z3;2;z3;z3;z3;z3];  %Dist-Body 1 to 2
SJDT(:,4)=[10;1;0;z3;z3;0;uy;ux;uy;ux];  %RotDr-Body 1 to Ground

%Initial generalized coordinate estimate
r10=[0;0;0];
p10=[1;0;0;0];
r20=[0;3;0];
p20=[1;0;0;0];
q0e=[r10;p10;r20;p20];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5   %Spatial Slider-Crank, Bod 1 RotDr
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nd=ngc-nc;

%SJDT(22,nh): Spatial Joint Data Table 
%SJTd(:,k)=[T;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; T=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph, 9=DistDr, 10=RotDr); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[4;1;0;z3;0.1*uy+0.12*uz;0;uz;ux;uz;ux];   %Rev-Body1 to Ground
SJDT(:,2)=[5;2;0;z3;z3;0;uz;ux;uz;ux];              %Tran-Body2 to Ground
SJDT(:,3)=[1;1;2;0.08*uz;z3;0.24;z3;z3;z3;z3];  %Dist-Body 1 to 2
SJDT(:,4)=[10;1;0;z3;0.1*uy+0.12*uz;0;uz;ux;uz;ux];  %RotDr-Body 1 to Ground

%Initial generalized coordinate estimate
r10=[0;0.1;0.12];
p10=[1;0;0;0];
r20=[0.25;0;0];
p20=[1;0;0;0];
q0e=[r10;p10;r20;p20];

end

end

