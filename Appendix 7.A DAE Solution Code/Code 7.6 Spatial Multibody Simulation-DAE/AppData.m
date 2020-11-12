function [nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q0,qd0]=...
    AppData(app)

if app==1   %Pendulum, Spherical to Ground
    
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
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;uxipr;uzipr;uxjpr;uzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; uxipr, uzipr, uxjpr, uzjpr
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
r10=[0;0;-1];
p10=[0;ux];
q0=[r10;p10];
omeg1pr0=[2;0;0];
r1d0=ATran(p10)*atil(omeg1pr0)*uz;
p1d0=0.5*GEval(p10)'*omeg1pr0;
qd0=[r1d0;p1d0];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2        %Top, Spherical to Ground
    
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
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;uxipr;uzipr;uxjpr;uzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; uxipr, uzipr, uxjpr, uzjpr
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
if app==3   %One Body Pendulum, Dist. to Ground
    
nb=1;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=1;       %Number of holonomic constraints
nhc=1;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
nv=ngc-nc;
nu=nc;
NTSDA=0;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;uxipr;uzipr;uxjpr;uzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; uxipr, uzipr, uxjpr, uzjpr
SJDT(:,1)=[1;1;0;zer;zer;1;zer;zer;zer;zer];   %Dist. - Body1 to Ground

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
r10=[1;0;0];
p10=[1;0;0;0];
q0=[r10;p10];
r1d0=[0;0;0];
p1d0=[0;0;0;0];
qd0=[r1d0;p1d0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==4   %Two Body Pendulum
    
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
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;uxipr;uzipr;uxjpr;uzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; uxipr, uzipr, uxjpr, uzjpr
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
if app==5   %One Body Cylindrical with Spring
    
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
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;uxipr;uzipr;uxjpr;uzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; uxipr, uzipr, uxjpr, uzjpr
SJDT(:,1)=[3;1;0;zer;zer;0;ux;uz;ux;uz];     %Cyl - Body1 to Ground

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
if app==6   %Spatial Slider-Crank
    
nb=2;       %Number of bodies
ngc=7*nb;   %number of generalized coordinates
nh=3;       %Number of holonomic constraints
nhc=11;      %Number of holonomic constraint equations
nc=nhc+nb;  %Number of constraint equations
NTSDA=0;    %Number of TSDA force elements

ux=[1;0;0];
uy=[0;1;0];
uz=[0;0;1];
zer=zeros(3,1);

%SJDT(22,nh): Joint Data Table 
%SJTd(:,k)=[t;i;j;sipr;sjpr;d;uxipr;uzipr;uxjpr;uzjpr]; 
    %k=joint No.; t=joint type(1=Dist,2=Sph,3=Cyl, 4=Rev, 5=Tran, 
    %6=Univ, 7=Strut, 8=Rev-Sph); i&j=bodies conn.,i>0;
    %si&jpr=vectors to Pi&j; d=dist.; uxipr, uzipr, uxjpr, uzjpr
SJDT(:,1)=[4;1;0;zer;0.1*ux+0.12*uy;0;ux;uz;ux;uz];   %Cyl-Body1 to Ground
SJDT(:,2)=[5;2;0;zer;zer;0;ux;uz;ux;uz];              %Tran-Body2 to Ground
SJDT(:,3)=[1;1;2;0.08*uy;0.02*uy;0.23;zer;zer;zer;zer];  %Dist-Body 1 to 2


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
r10=[0.1;0.12;0];
p10=[1;0;0;0];
r20=[0;0;0.1027];
p20=[1;0;0;0];
q0=[r10;p10;r20;p20];
r1d0=[0;0;0];
p1d0=0.5*EEval(p10)'*120*uz;
r2d0=[0;0;9.378];
p2d0=[0;0;0;0];
qd0=[r1d0;p1d0;r2d0;p2d0];

end

end

