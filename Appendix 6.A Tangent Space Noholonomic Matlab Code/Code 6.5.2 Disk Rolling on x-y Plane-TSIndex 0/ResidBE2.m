function R=ResidBE2(t,yk,k,Lam,ue,B,H,q0,qd0,...
    V,U,W,X,par,utol,Btol,Htol)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

y=yk+0.5*h*k;

v=[y(1);y(2);y(3);y(4);y(5)];
w=[y(6);y(7);y(8)];

[u,uiter]=usolv(t,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Bcorr(t,q,B,U,par);
C=Ceval(t,q,par);
[H,Hiter]= Hcorr(H,X,C,par);
I9=eye(9);
D=(I9-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(t,q,par);
qd=D*w+(I9-X*H*C)*qd0+X*H*Nu;
M=Meval(q,par);
Gam=Gameval(t,q,qd,par);
S=Seval(q,qd,par);
QA=QAEval(t,q,qd,par);

kv=[k(1);k(2);k(3);k(4);k(5)];
kw=[k(6);k(7);k(8)];

R1=kv-V'*(D*w+(I9-X*H*C)*qd0+X*H*Nu);
R2=M*D*kw+C'*Lam-M*X*H*Gam-S-QA;

R=[R1;R2];


end

