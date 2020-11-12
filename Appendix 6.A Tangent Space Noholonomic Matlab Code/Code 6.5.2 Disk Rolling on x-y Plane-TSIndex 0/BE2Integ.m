function[v,vd,w,wd,Lam,jiter,RM,k1,k2,RR]=...
    BE2Integ(tnm,n,nRepar,tn,Vv,Vvd,Ww,Wwd,LLam,Uu,Kk1,Kk2,q0,qd0,...
    V,U,W,X,B,H,h,hmax,utol,intol,Btol,Htol,par,J2,errcontr)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

 %First Order BE2 Coefficients

a=1/2;
A=0.5*[1,0;1,1];
c=[1/2;1];
b=0.5*[1,1];

% Start Integration Step
Lam=LLam(:,n-1);
vnm=Vv(:,n-1);
vdnm=Vvd(:,n-1);
wnm=Ww(:,n-1);
wdnm=Wwd(:,n-1);

%Jacobian Evaluation
ue=Uu(:,n-1);
[u,uiter]=usolv(tnm,ue,vnm,q0,V,U,B,par);
q=q0+V*vnm-U*u;
M=Meval(q,par);
C=Ceval(tnm,q,par);
[H,Hiter]= Hcorr(H,X,C,par);
I9=eye(9);
D=(I9-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(tnm,q,par);
qd=D*wnm+(I9-X*H*C)*qd0+X*H*Nu;

J2=J2Eval(tnm,q,qd,qd0,wnm,wdnm,Lam,U,V,W,X,B,H,par);

Jac=[eye(nv),zeros(nv,nw),zeros(nv,nh+nd);zeros(nq,nv),M*D,C']+(h/2)*J2;
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

t1=tnm+(h/2);
y1k=[vnm;wnm];

while err>intol 
    
y1=y1k+h*a*k1;

R1=ResidBE2(t1,y1k,k1,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);

if i==1
    NR1=norm(R1);
    R10=R1;
end

% Newton Correction

z=-Jac\R1;
delk1=[z(1);z(2);z(3);z(4);z(5);z(6);z(7);z(8)];
delLam=[z(9);z(10);z(11);z(12);z(13);z(14)];
k1=k1+delk1;
Lam=Lam+delLam;

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

t2=tnm+(h);
y2k=y1k+h*A(2,1)*k1;

while err>intol
    
y2=y2k+h*a*k2;
% Residual Calculation

R2=ResidBE2(t2,y2k,k2,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol);
if i==1
    NR2=norm(R2);
    R20=R2;
end

% Newton Correction

z=-Jac\R2;
delk2=[z(1);z(2);z(3);z(4);z(5);z(6);z(7);z(8)];
delLam=[z(9);z(10);z(11);z(12);z(13);z(14)];
k2=k2+delk2;
Lam=Lam+delLam;

err=norm(R2);

i=i+1;
end
j2=i-1;

R=[NR1;NR2];
RR=[R10,R20];
jiter=j1;
kk=[k1,k2];


%Evaluate y, v, w,vd, wd

y=y1k+h*kk*b';
v=[y(1);y(2);y(3);y(4);y(5)];
w=[y(6);y(7);y(8)];

tn=tnm+h;
[u,uiter]=usolv(tn,u,v,q0,V,U,B,par);
[f,Lam,RM]=ODEfunct(tn,y,u,V,U,W,X,B,H,q0,qd0,par);
vd=[f(1);f(2);f(3);f(4);f(5)];
wd=[f(6);f(7);f(8)];


end

