function [f,Lam,RM]=ODEfunct(t,y,ue,V,U,W,X,B,H,q0,qd0,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ,h,utol,Btol,Htol]...
    =Partpar(par);

v=[y(1);y(2);y(3);y(4)];
w=y(5);
[u,uiter]=usolv(t,ue,v,q0,V,U,B,utol,par);
q=q0+V*v-U*u;
Phiq=Phiqeval(t,q,par);
[B,Biter]=Bcorr(t,q,B,U,Btol,par);
C=Ceval(t,q,par);
[H,Hiter]= Hcorr(H,X,C,Htol,par);
D=(eye(nq)-X*H*C)*W;
[Nu,Nusq]=NuNusqeval(t,q,par);
qd=D*w+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Gam=Gameval(t,q,qd,par);
M=Meval(q,par);
QA=QAEval(t,q,qd,par);
S=Seval(q,qd,par);

RM=[eye(4),zeros(4,1),zeros(4,5);zeros(6,4),M*D,C'];
RHS=[V'*(D*w+(eye(nq)-X*H*C)*qd0+X*H*Nu);M*X*H*Gam+S+QA];

z=RM\RHS;

vd=[z(1);z(2);z(3);z(4)];
wd=z(5);
f=[vd;wd];
Lam=[z(6);z(7);z(8);z(9);z(10)];



end



