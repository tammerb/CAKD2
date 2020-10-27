function [nb,ngc,nh,nc,nv,nu,NTSDA,SJDT,SMDT,STSDAT,q0,qd0]=...
    AppData(app)

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
z3=zeros(3,1);
z4=zeros(4,1);

if app==1   %Body in Cylindrical with Ground
    
nb=1;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates

NTSDA=1;    %Number of TSDA force elements

%SJDT(28,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
    %a,b or R=joint dimensions; mus,mud=coeffidients of stat. &dyn. frict.
    %ms=address of firsrt Lag. Mult,;nm= No. of Constr. Eq.
SJDT(:,1)=[3;1;0;z3;z3;0;ux;uz;ux;uz;0.1;0.05;0.24;0.2;1;4];%Cyl - Body1 to Ground

nh=1;       %Number of holonomic constraints
nhc=4;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[1;0.1;0.1;0.1]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force
STSDAT(:,1)=[1;0;uy+uz;-ux+uy+uz;10;0;1;0];

%Initial generalized coordinates
r10=[0;0;0];
p10=[1;0;0;0];
q0=[r10;p10];
r1d0=[0;0;0];
p1d0=0.5*EEval(p10)'*10*uz;
qd0=[r1d0;p1d0];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Spatial Slider-Crank
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
NTSDA=0;    %Number of TSDA force elements

%SJDT(28,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
    %a,b or R=joint dimensions; mus,mud=coeffidients of stat. &dyn. frict.
    %ms=address of firsrt Lag. Mult,;nm= No. of Constr. Eq.
SJDT(:,1)=[4;1;0;z3;0.1*uy+0.12*uz;0;uz;ux;uz;ux;0.1;0.1;0.1;0.08;1;5];   %Rev-Body1 to Ground
SJDT(:,2)=[5;2;0;z3;z3;0;uz;ux;uz;ux;0.1;0.1;0;0;6;5];              %Tran-Body2 to Ground
SJDT(:,3)=[1;1;2;0.08*uz;0.02*uz;0.3;z3;z3;z3;z3;0.1;0.1;0.1;0.08;11;1];%Dist-Body 1 to 2

nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==3   %Four Body Slider
    
nb=4;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
NTSDA=4;    %Number of TSDA force elements

%SJDT(28,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
    %a,b or R=joint dimensions; mus,mud=coeffidients of stat. &dyn. frict.
    %ms=address of firsrt Lag. Mult,;nm= No. of Constr. Eq.
SJDT(:,1)=[5;1;0;z3;z3;0;uz;ux;uz;ux;0.1;0.1;0.12;0.1;1;5];   %Tran. - Body1 to Ground
SJDT(:,2)=[5;2;0;z3;z3;0;ux;uy;ux;uy;0.1;0.1;0.12;0.1;6;5];   %Tran. - Body2 to Ground
SJDT(:,3)=[5;3;0;z3;z3;0;ux;uz;ux;uz;0.1;0.1;0.12;0.1;11;5];   %Tran. - Body3 to Ground
SJDT(:,4)=[5;4;0;z3;z3;0;uz;ux;uz;ux;0.1;0.1;0.12;0.1;16;5];   %Tran. - Body3 to Ground
SJDT(:,5)=[1;1;2;z3;z3;5;z3;z3;z3;z3;0;0;0;0;21;1];   %Dist. - Body1 to 2
SJDT(:,6)=[1;2;3;z3;z3;7;z3;z3;z3;z3;0;0;0;0;22;1];   %Dist. - Body1 to 2

nh=6;       %Number of holonomic constraints
nhc=22;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[2;1;1;1],[2;1;1;1],[6;1;1;1],[6;1;1;1]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force
STSDAT(:,1)=[1;0;z3;-10*ux;2;0;13;0];
STSDAT(:,2)=[2;0;z3;-10*uy;2;0;14;0];
STSDAT(:,3)=[3;0;z3;-10*uz;2;0;13;0];
STSDAT(:,4)=[1;4;z3;10*ux;10;0;11;0];
%Initial generalized coordinates
r10=[4;0;0];
p10=[1;0;0;0];
q10=[r10;p10];
r1d0=[0;0;0];
p1d0=[0;0;0;0];
q1d0=[r1d0;p1d0];
r20=[0;3;0];
p20=[1;0;0;0];
q20=[r20;p20];
r2d0=[0;0;0];
p2d0=[0;0;0;0];
q2d0=[r2d0;p2d0];
r30=[0;0;6.32];
p30=[1;0;0;0];
q30=[r30;p30];
r3d0=[0;0;0];
p3d0=[0;0;0;0];
q3d0=[r3d0;p3d0];
r40=[5;0;0];
p40=[1;0;0;0];
q40=[r40;p40];
r4d0=[-1;0;0];
p4d0=[0;0;0;0];
q4d0=[r4d0;p4d0];
q0=[q10;q20;q30;q40];
qd0=[q1d0;q2d0;q3d0;q4d0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Rotating Disk with Translating Body
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates

NTSDA=1;    %Number of TSDA force elements

%SJDT(28,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
    %a,b or R=joint dimensions; mus,mud=coeffidients of stat. &dyn. frict.
    %ms=address of firsrt Lag. Mult,;nm= No. of Constr. Eq.
SJDT(:,1)=[4;1;0;z3;z3;0;ux;uz;ux;uz;0.1;0.1;0.08;0.07;1;5];%Rev - 1 to Ground
SJDT(:,2)=[5;1;2;ux;z3;0;ux;uy;ux;uy;0.1;0.1;0.08;0.07;6;5];%Tran. - 1 to 2

nh=2;       %Number of holonomic constraints
nhc=10;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;

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
STSDAT(:,1)=[1;2;ux+10*uy;z3;10;0;10.1;0];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==5   %One Body Pendulum, Cyl about x axis
    
nb=1;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
NTSDA=0;    %Number of TSDA force elements

%SJDT(28,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
    %a,b or R=joint dimensions; mus,mud=coeffidients of stat. &dyn. frict.
    %ms=address of firsrt Lag. Mult,;nm= No. of Constr. Eq.
    SJDT(:,1)=[3;1;0;-uz;z3;0;uz;ux;uz;ux;0.1;0.1;0.3;0.25;1;4];  %Rev. Bod1-Grnd
   

nh=1;       %Number of holonomic constraints
nhc=4;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[75;30;30;30]];

%STSDAT(12,1): TSDA Data Table
if NTSDA==0
STSDAT=zeros(12,NTSDA);
end
%STSDAT(:,T)=[i;j;sipr;sjpr;K;C;el0;F];  
    %T=TSDA No.; i&j=bodies conn.;si&jpr=vectors to Pi&j; K=spring constant;
    %C=damping coefficient; el0=spring free length; F=const. force

%Initial generalized coordinates
r10=[0;0;-1];
p10=[1;0;0;0];
q0=[r10;p10];
r1d0=[0.1;0;0];
p1d0=0.5*EEval(p10)'*1*ux;
qd0=[r1d0;p1d0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==6   %One Body Translation Along x axis
    
nb=1;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
NTSDA=0;    %Number of TSDA force elements

%SJDT(28,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;vxipr;vzipr;vxjpr;vzjpr;a;b or R;mus;mud;mc;nm]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; vxipr, vzipr, vxjpr, vzjpr;
    %a,b or R=joint dimensions; mus,mud=coeffidients of stat. &dyn. frict.
    %ms=address of firsrt Lag. Mult,;nm= No. of Constr. Eq.
    SJDT(:,1)=[5;1;0;z3;z3;0;uz;ux;uz;ux;0.1;0.1;0.3;0.25;1;4];  %Tran. Bod1-Grnd
   

nh=1;       %Number of holonomic constraints
nhc=5;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;

%SMDT(4,nb): Mass Data Table (With diagonal inertia matrix) 
%SMDT=[[m1;J11;J12,J13],...,[mnb;Jnb1;Jnb2;Jnb3]]
SMDT=[[75;30;30;30]];

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
q0=[r10;p10];
r1d0=[1;0;0];
p1d0=z4;
qd0=[r1d0;p1d0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

