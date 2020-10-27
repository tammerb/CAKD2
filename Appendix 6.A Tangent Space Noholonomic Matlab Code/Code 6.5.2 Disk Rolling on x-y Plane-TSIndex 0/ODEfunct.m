function [f,Lam,RM]=ODEfunct(t,y,ue,V,U,W,X,B,H,q0,qd0,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

v=[y(1);y(2);y(3);y(4);y(5)];
w=[y(6);y(7);y(8)];

[u,uiter]=usolv(t,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
Phiq=Phiqeval(t,q,par);
[B,Biter]=Bcorr(t,q,B,U,par);
C=Ceval(t,q,par);
[H,Hiter]= Hcorr(H,X,C,par);
D=(eye(nq)-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(t,q,par);
qd=D*w+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Gam=Gameval(t,q,qd,par);
M=Meval(q,par);
QA=QAEval(t,q,qd,par);
S=Seval(q,qd,par);

RM=[eye(nv),zeros(nv,nw),zeros(nv,nh+nd);zeros(nq,nv),M*D,C'];
RHS=[V'*(D*w+(eye(nq)-X*H*C)*qd0+X*H*Nu);M*X*H*Gam+S+QA];

z=RM\RHS;

vd=[z(1);z(2);z(3);z(4);z(5)];
wd=[z(6);z(7);z(8)];
f=[vd;wd];
Lam=[z(9);z(10);z(11);z(12);z(13);z(14)];



end



