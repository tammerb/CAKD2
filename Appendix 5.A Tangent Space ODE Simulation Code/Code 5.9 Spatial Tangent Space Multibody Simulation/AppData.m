function [nb,ngc,nh,nc,nv,nu,NTSDA,SJDT,SMDT,STSDAT,q0,qd0]=...
    AppData(app)

z3=[0;0;0];
z4=[0;0;0;0];
I3=eye(3);
ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];

if app==1        %Spin Stabilized Top
    
nb=1;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=1;       %Number of holonomic constraints
nhc=3;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[2;1;0;-uz;zer;0;zer;zer;zer;zer];    %Sph Jt - Body1 and ground

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix)
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[30;90;90;30];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force

%Initial generalized coordinates
r10=[0;0;1];
p10=[1;0;0;0];
q0=[r10;p10];
omeg1pr0=[10^-12;10^-12;13.5];
r1d0=atil(omeg1pr0)*uz;
p1d0=0.5*GEval(p10)'*omeg1pr0;
qd0=[r1d0;p1d0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Spatial Double Pendulum
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=2;       %Number of holonomic constraints
nhc=4;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
    SJDT(:,1)=[2;1;0;zer;zer;0;zer;zer;zer;zer];  %Sph. - Body 1 to Ground
    SJDT(:,2)=[1;1;2;-uz;uz;1;zer;zer;zer;zer];   %Dist. - Body 1 to 2

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[75;30;30;30],[75;30;30;30]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force

%Initial generalized coordinates
r10=[0;0;0];
p10=[1;0;0;0];
r20=[0;0;-3];
p20=[1;0;0;0];
q0=[r10;p10;r20;p20];
r1d0=[0;0;0];
p1d0=[0;0;0;0];
r2d0=[0;0;0];
p2d0=0.5*GEval(p20)'*[0;5;0];
qd0=[r1d0;p1d0;r2d0;p2d0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %One Body in Cylindrical with Spring
    
nb=1;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=1;       %Number of holonomic constraints
nhc=4;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=1;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[3;1;0;zer;zer;0;ux;uz;ux;uz];     %Cyl - Body1 to Ground

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[1;1;1;1]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force
STSDAT(:,1)=[1;0;uy+uz;-ux+uy+uz;100;0;1;0];

%Initial generalized coordinates
r10=[0;0;0];
p10=[1;0;0;0];
q0=[r10;p10];
r1d0=[0;0;0];
p1d0=0.5*EEval(p10)'*10*uz;
qd0=[r1d0;p1d0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Spatial Slider-Crank
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[4;1;0;zer;0.1*uy+0.12*uz;0;uz;ux;uz;ux];   %Rev-Body1 to Ground
SJDT(:,2)=[5;2;0;zer;zer;0;uz;ux;uz;ux];              %Tran-Body2 to Ground
SJDT(:,3)=[1;1;2;0.08*uz;0.02*uz;0.3;zer;zer;zer;zer];  %Dist-Body 1 to 2

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[0.5;0.2;0.2;0.2],[5;0.2;0.2;0.2]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force

%Initial generalized coordinates
r10=[0;0.1;0.12];
p10=[1;0;0;0];
r20=[0.2182;0;0];
p20=[1;0;0;0];
q0=[r10;p10;r20;p20];
r1d0=[0;0;0];
p1d0=0.5*EEval(p10)'*120*ux;
r2d0=[4.406;0;0];
p2d0=[0;0;0;0];
qd0=[r1d0;p1d0;r2d0;p2d0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5   %4-Translating Mass Model
    
nb=4;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=6;       %Number of holonomic constraints
nhc=22;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=4;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[5;1;0;zer;zer;0;uz;ux;uz;ux];   %Tran-Body1 to Ground
SJDT(:,2)=[5;2;0;zer;zer;0;uz;uy;uz;uy];   %Tran-Body2 to Ground
SJDT(:,3)=[5;3;0;zer;zer;0;ux;uz;ux;uz];   %Tran-Body3 to Ground
SJDT(:,4)=[5;4;0;zer;zer;0;uz;ux;uz;ux];   %Tran-Body4 to Ground
SJDT(:,5)=[1;1;2;zer;zer;5;zer;zer;zer;zer];  %Dist-Body 1 to 2
SJDT(:,6)=[1;2;3;zer;zer;7;zer;zer;zer;zer];  %Dist-Body 2 to 3

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[2;0.2;0.2;0.2],[2;0.2;0.2;0.2],[6;0.2;0.2;0.2],[6;0.2;0.2;0.2]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force
STSDAT(:,1)=[1;0;10*ux;zer;2;0;14;0];       %mass1 to ground
STSDAT(:,2)=[2;0;10*uy;zer;2;0;13;0];       %mass2 to ground
STSDAT(:,3)=[3;0;10*uz;zer;2;0;16.32;0];       %mass3 to ground
STSDAT(:,4)=[1;4;zer;10*ux;10;0;11;0];       %mass4 to mass1

%Initial generalized coordinates
r10=[4;0;0];
p10=[1;0;0;0];
r20=[0;3;0];
p20=[1;0;0;0];
r30=[0;0;6.32];
p30=[1;0;0;0];
r40=[5;0;0];
p40=[1;0;0;0];
q0=[r10;p10;r20;p20;r30;p30;r40;p40];
r1d0=[0;0;0];
p1d0=[0;0;0;0];
r2d0=[0;0;0];
p2d0=[0;0;0;0];
r3d0=[0;0;0];
p3d0=[0;0;0;0];
r4d0=[-1;0;0];
p4d0=[0;0;0;0];
qd0=[r1d0;p1d0;r2d0;p2d0;r3d0;p3d0;r4d0;p4d0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==6   %Rotating Disk with Translating Body
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=2;       %Number of holonomic constraints
nhc=10;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=1;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr
SJDT(:,1)=[4;1;0;z3;z3;0;ux;uz;ux;uz];  %Rev. - Body 1 to Ground
SJDT(:,2)=[5;1;2;ux;z3;0;ux;uy;ux;uy];   %Tran. - Body 1 to 2

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[10;10;10;10],[5;5;5;5]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force
STSDAT(:,1)=[1;2;ux+10*uy;z3;10;1;10.1;0];

%Initial generalized coordinates
r10=[0;0;0];
p10=[1;0;0;0];
r20=[1;0;0];
p20=[1;0;0;0];
q0=[r10;p10;r20;p20];
r1d0=[0;0;0];
p1d0=[0;0;0;0];
r2d0=[0;0;0];
p2d0=[0;0;0;0];
qd0=[r1d0;p1d0;r2d0;p2d0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

