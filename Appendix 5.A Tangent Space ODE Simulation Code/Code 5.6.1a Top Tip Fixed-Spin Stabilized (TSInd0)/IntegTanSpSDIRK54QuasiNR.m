function[v,vd,vdd,Lam,jiter,R,k1,k2,k3,k4,k5,RR,h]=IntegTanSpSDIRK54QuasiNR...
    (n,npar,tn,Vv,Vvd,Vvdd,Qdd,LLam,Uu,Kk1,Kk2,Kk3,Kk4,Kk5,q0,V,U,B,...
    h,h2,utol,intol,Btol,par,Jinv)

 %Runge-Kutta Coefficients

a=1/4;
aa=0.0625;
As=[1/4,0,0,0,0;1/2,1/4,0,0,0;17/50,-1/25,1/4,0,0;...
    371/1360,-137/2720,15/544,1/4,0;...
    25/24,-49/48,125/16,-85/12,1/4];
c=[1/4;3/4;11/20;1/2;1];
b=[25/24,-49/48,125/16,-85/12,1/4];
d=[59/48,-17/96,225/32,-85/12,0];
AA=As*As;
bbar=b*As;
D=d*As;

% Start Integration Step
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
vddnm=Vvdd(:,n-1);
qddnm=Qdd(:,n-1);
u=Uu(:,n-1);

% Solve Discretized Equations

% Stage 1
i=1;        %Set solution iteration counter
err=intol+1;

if n==2
    qdd=qddnm;
    Lam=LLam(:,n-1);
end

if n>2
    qdd=2*qddnm-Qdd(:,n-2);
    Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

vdd=(V')*qdd;

if n-npar==1
k1=vdd;
Lam=LLam(:,n-1);
end
if n-npar>1
k1=2*Kk1(:,n-1)-Kk1(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

z1k=vnm+h*c(1)*vdnm;
z1dk=vdnm;
v=z1k+h2*aa*k1;
vd=z1dk+h*a*k1;

while err>intol;    

R1=Resid(v,vd,k1,Lam,u,par,q0,V,U,B,utol,Btol,h);

if i==1;
    NR1=norm(R1);
    R10=R1;
end
% Quasi Newton Correction
z=-Jinv*R1;
delk1=[z(1);z(2);z(3)];
delLam=[z(4);z(5);z(6);z(7)];
k1=k1+delk1;
Lam=Lam+delLam;
v=z1k+h2*aa*k1;
vd=z1dk+h*a*k1;
err=norm(R1);
i=i+1;
end
j1=i-1;

% Stage 2
i=1;        %Set solution iteration counter
err=intol+1;

if n-npar==1
k2=k1;
Lam=LLam(:,n-1);
end
if n-npar>1
k2=2*Kk2(:,n-1)-Kk2(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

z2k=vnm+h*c(2)*vdnm+h2*AA(2,1)*k1;
z2dk=vdnm+h*As(2,1)*k1;
v=z2k+h2*aa*k2;
vd=z2dk+h*a*k2;
while err>intol;    
% Residual Calculation

R2=Resid(v,vd,k2,Lam,u,par,q0,V,U,B,utol,Btol,h);
if i==1;
    NR2=norm(R2);
    R20=R2;
end
% Quasi Newton Correction
z=-Jinv*R2;
delk2=[z(1);z(2);z(3)];
delLam=[z(4);z(5);z(6);z(7)];
k2=k2+delk2;
Lam=Lam+delLam;
v=z2k+h2*aa*k2;
vd=z2dk+h*a*k2;
err=norm(R2);

i=i+1;
end
j2=i-1;

% Stage 3
i=1;        %Set solution iteration counter
err=intol+1;

if n-npar==1
k3=k2;
Lam=LLam(:,n-1);
end
if n-npar>1
k3=2*Kk3(:,n-1)-Kk3(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

z3k=vnm+h*c(3)*vdnm+h2*(AA(3,1)*k1+AA(3,2)*k2);
z3dk=vdnm+h*(As(3,1)*k1+As(3,2)*k2);
v=z3k+h2*aa*k3;
vd=z3dk+h*a*k3;
while err>intol;    

R3=Resid(v,vd,k3,Lam,u,par,q0,V,U,B,utol,Btol,h);
if i==1;
    NR3=norm(R3);
    R30=R3;
end
% Quasi Newton Correction
z=-Jinv*R3;
delk3=[z(1);z(2);z(3)];
delLam=[z(4);z(5);z(6);z(7)];
k3=k3+delk3;
Lam=Lam+delLam;
v=z3k+h2*aa*k3;
vd=z3dk+h*a*k3;
err=norm(R3);
i=i+1;
end
j3=i-1;

% Stage 4
i=1;        %Set solution iteration counter
err=intol+1;

if n-npar==1
k4=k3;
Lam=LLam(:,n-1);
end
if n-npar>1
k4=2*Kk4(:,n-1)-Kk4(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

z4k=vnm+h*c(4)*vdnm+h2*(AA(4,1)*k1+AA(4,2)*k2+AA(4,3)*k3);
z4dk=vdnm+h*(As(4,1)*k1+As(4,2)*k2+As(4,3)*k3);
v=z4k+h2*aa*k4;
vd=z4dk+h*a*k4;
while err>intol;    

R4=Resid(v,vd,k4,Lam,u,par,q0,V,U,B,utol,Btol,h);
if i==1;
    NR4=norm(R4);
    R40=R4;
end
% Quasi Newton Correction
z=-Jinv*R4;
delk4=[z(1);z(2);z(3)];
delLam=[z(4);z(5);z(6);z(7)];
k4=k4+delk4;
Lam=Lam+delLam;
v=z4k+h2*aa*k4;
vd=z4dk+h*a*k4;
err=norm(R4);
%err=sqrt(R'*R)
i=i+1;
end
j4=i-1;

% Stage 5
i=1;        %Set solution iteration counter
err=intol+1;

if n-npar==1
k5=k4;
Lam=LLam(:,n-1);
end
if n-npar>1
k5=2*Kk5(:,n-1)-Kk5(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

z5k=vnm+h*c(5)*vdnm+...
    h2*(AA(5,1)*k1+AA(5,2)*k2+AA(5,3)*k3+AA(5,4)*k4);
z5dk=vdnm+h*(As(5,1)*k1+As(5,2)*k2+As(5,3)*k3+As(5,4)*k4);
v=z5k+h2*aa*k5;
vd=z5dk+h*a*k5;
while err>intol;    
    
R5=Resid(v,vd,k5,Lam,u,par,q0,V,U,B,utol,Btol,h);
if i==1;
    NR5=norm(R5);
    R50=R5;
end
% Quasi Newton Correction
z=-Jinv*R5;
delk5=[z(1);z(2);z(3)];
delLam=[z(4);z(5);z(6);z(7)];
k5=k5+delk5;
Lam=Lam+delLam;
v=z5k+h2*aa*k5;
vd=z5dk+h*a*k5;
err=norm(R5);
i=i+1;
end
j5=i-1;

R=[NR1;NR2;NR3;NR4;NR5];
RR=[R10,R20,R30,R40,R50];
jiter5=[j1;j2;j3;j4;j5];
jiter=j1;
kk=[k1,k2,k3,k4,k5];


%Evaluate v and vd and vdd
v=vnm+h*vdnm+h2*kk*bbar';
vd=vdnm+h*kk*b';

[u,Iteru]=usolv(u,v,q0,V,U,B,utol);
[f,Lam,B,E]=ODEfunct(tn,v,vd,u,V,U,B,q0,utol,Btol,par);
vdd=f;



end