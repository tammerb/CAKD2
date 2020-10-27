function [f,RM]=ODEfunct(t,y,ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,L,Jpr)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
[v,w] = yPart(y);
[u,uiter] = usolv(t,ue,v,q0,V,U,B,utol,par,L);
q=q0+V*v-U*u;
Phiq=Phiqeval(t,q,par,L);
[B,Biter]=Bcorr(t,q,B,U,Btol,par,L);
C=Ceval(t,q,par,L);
[H,Hiter]=Hcorr(H,X,C,Htol,par);
D2=(eye(nq)-X*H*C)*W;
Nu=Nueval(t,q,par);
qd=D2*w+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Gam=Gameval(t,q,qd,par,L);
M=Meval(q,par,Jpr);
RM=D2'*M*D2;
QA=QAEval(t,q,qd,par);
S=Seval(q,qd,par,Jpr);
RHS=D2'*(M*X*H*Gam+S+QA);
wd=RM\RHS;
Nu=Nueval(t,q,par);
vd=V'*D2*w+V'*(eye(nq)-X*H*C)*qd0+V'*X*H*Nu;
f=[vd;wd];


end



