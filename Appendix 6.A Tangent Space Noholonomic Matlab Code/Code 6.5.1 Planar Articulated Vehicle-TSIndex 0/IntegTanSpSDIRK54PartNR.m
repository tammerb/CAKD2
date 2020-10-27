function[v,vd,w,wd,Lam,jiter,RM,k1,k2,k3,k4,k5,RR,h]=...
    IntegTanSpSDIRK54PartNR(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,Kk2,...
    Kk3,Kk4,Kk5,q0,qd0,V,U,W,X,B,H,h,hmax,utol,intol,Btol,Htol,par,...
    J2,errcontr)

[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ,h,utol,Btol,Htol]...
    =Partpar(par);

 %First Order SDIRK Coefficients

a=1/4;
A=[1/4,0,0,0,0;1/2,1/4,0,0,0;17/50,-1/25,1/4,0,0;...
    371/1360,-137/2720,15/544,1/4,0;...
    25/24,-49/48,125/16,-85/12,1/4];
c=[1/4;3/4;11/20;1/2;1];
b=[25/24,-49/48,125/16,-85/12,1/4];
d=[59/48,-17/96,225/32,-85/12,0];

% Start Integration Step
Lam=LLam(:,n-1);
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
wnm=Ww(:,n-1);
wdnm=Wwd(:,n-1);

%Jacobian Evaluation
ue=Uu(:,n-1);
[u,uiter]=usolv(tnm,ue,vnm,q0,V,U,B,utol,par);
q=q0+V*vnm-U*u;
M=Meval(q,par);
C=Ceval(tnm,q,par);
[H,Hiter]= Hcorr(H,X,C,Htol,par);
I6=eye(6);
D=(I6-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(tnm,q,par);
qd=D*wnm+(I6-X*H*C)*qd0+X*H*Nu;

J2=J2Eval(tnm,q,qd,wnm,wdnm,Lam,U,V,W,X,B,H,q0,qd0,par);

Jac=[eye(4),zeros(4,1),zeros(4,5);zeros(6,4),M*D,C']+(h/4)*J2;
RM=Jac;

% Solve Discretized Equations

% Stage 1
i=1;        %Set solution iteration counter
err=intol+1;

if n-nRepar==1
k1=[vdnm;wdnm];
Lam=LLam(:,n-1);
end
if n-nRepar>1
k1=2*Kk1(:,n-1)-Kk1(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

t1=tnm+(h/4);
y1k=[vnm;wnm];
y1=y1k+h*a*k1;

while err>intol;    

R1=ResidSDK54(t1,y1k,k1,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);

if i==1;
    NR1=norm(R1);
    R10=R1;
end

% Newton Correction

z=-Jac\R1;
delk1=[z(1);z(2);z(3);z(4);z(5)];
delLam=[z(6);z(7);z(8);z(9);z(10)];
k1=k1+delk1;
Lam=Lam+delLam;
y1=y1k+h*a*k1;
err=norm(R1);
i=i+1;
end
j1=i-1;

% Stage 2
i=1;        %Set solution iteration counter
err=intol+1;

if n-nRepar==1
k2=k1;
Lam=LLam(:,n-1);
end
if n-nRepar>1
k2=2*Kk2(:,n-1)-Kk2(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

t2=tnm+(3*h/4);
y2k=y1k+h*A(2,1)*k1;
y2=y2k+h*a*k2;
while err>intol;    
% Residual Calculation

R2=ResidSDK54(t2,y2k,k2,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);
if i==1;
    NR2=norm(R2);
    R20=R2;
end

% Newton Correction

z=-Jac\R2;
delk5=[z(1);z(2);z(3);z(4);z(5)];
delLam=[z(6);z(7);z(8);z(9);z(10)];
k2=k2+delk5;
Lam=Lam+delLam;
y2=y2k+h*a*k2;
err=norm(R2);

i=i+1;
end
j2=i-1;

% Stage 3
i=1;        %Set solution iteration counter
err=intol+1;

if n-nRepar==1
k3=k2;
Lam=LLam(:,n-1);
end
if n-nRepar>1
k3=2*Kk3(:,n-1)-Kk3(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

t3=tnm+(11*h/20);
y3k=y1k+h*(A(3,1)*k1+A(3,2)*k2);
y3=y3k+h*a*k3;
while err>intol;    

R3=ResidSDK54(t3,y3k,k3,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);
if i==1;
    NR3=norm(R3);
    R30=R3;
end

% Newton Correction

z=-Jac\R3;
delk3=[z(1);z(2);z(3);z(4);z(5)];
delLam=[z(6);z(7);z(8);z(9);z(10)];
k3=k3+delk3;
Lam=Lam+delLam;
y3=y3k+h*a*k3;
err=norm(R3);

i=i+1;
end
j3=i-1;

% Stage 4
i=1;        %Set solution iteration counter
err=intol+1;

if n-nRepar==1
k4=k3;
Lam=LLam(:,n-1);
end
if n-nRepar>1
k4=2*Kk4(:,n-1)-Kk4(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

t4=tnm+h/2;
y4k=y1k+h*(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3);
y4=y4k+h*a*k4;
while err>intol;    

R4=ResidSDK54(t4,y4k,k4,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);
if i==1;
    NR4=norm(R4);
    R40=R4;
end

% Newton Correction

z=-Jac\R4;
delk4=[z(1);z(2);z(3);z(4);z(5)];
delLam=[z(6);z(7);z(8);z(9);z(10)];
k4=k4+delk4;
Lam=Lam+delLam;
y4=y4k+h*a*k4;
err=norm(R4);

i=i+1;
end
j4=i-1;

% Stage 5
i=1;        %Set solution iteration counter
err=intol+1;

if n-nRepar==1
k5=k4;
Lam=LLam(:,n-1);
end
if n-nRepar>1
k5=2*Kk5(:,n-1)-Kk5(:,n-2);
Lam=2*LLam(:,n-1)-LLam(:,n-2);
end

t5=tnm+h;
y5k=y1k+h*(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3+A(5,4)*k4);
y5=y5k+h*a*k5;
while err>intol;    
    
R5=ResidSDK54(t5,y5k,k5,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);
if i==1;
    NR5=norm(R5);
    R50=R5;
end

% Newton Correction

z=-Jac\R5;
delk5=[z(1);z(2);z(3);z(4);z(5)];
delLam=[z(6);z(7);z(8);z(9);z(10)];
k5=k5+delk5;
Lam=Lam+delLam;
y5=y5k+h*a*k5;
err=norm(R5);


i=i+1;
end
j5=i-1;

R=[NR1;NR2;NR3;NR4;NR5];
RR=[R10,R20,R30,R40,R50];
jiter5=[j1;j2;j3;j4;j5];
jiter=j1;
kk=[k1,k2,k3,k4,k5];


%Evaluate y, v, w,vd, wd

y=y1k+h*kk*b';
v=[y(1);y(2);y(3);y(4)];
w=y(5);

tn=tnm+h;
[u,uiter]=usolv(tn,u,v,q0,V,U,B,utol,par);
[f,Lam,RM]=ODEfunct(tn,y,u,V,U,W,X,B,H,q0,qd0,par);
vd=[f(1);f(2);f(3);f(4)];
wd=f(5);

hopt=1;
if errcontr==1
% Error Control
%Evaluate yalt and error with Atol=0.00001, sc=Atol, and p=5

yalt=y1k+h*kk*d';

%Change step size
err=sqrt(((y-yalt)'*(y-yalt))/(0.0000000001*(nv+nw)));
hopt=h*(1/err)^(1/6);
h=0.5*hopt;

if h>hmax
    h=hmax; 
end
end

end