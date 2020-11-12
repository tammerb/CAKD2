function [f,RM]=ODEfunct(t,y,ue,V,U,W,X,B,H,q0,qd0,...
    utol,Btol,Htol,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[v,w] = yPart(y)
[u,uiter]=usolv(t,ue,v,q0,V,U,B,utol,par,J,dpP0,dpP1,dpP2,...
    bp,ux,uy,uz,P,A1,apppsa,atpppsa);
q=q0+V*v-U*u;
Phiq=P1(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
[B,Biter]=Bcorr(t,q,B,U,Btol,par,J,dpP0,dpP1,dpP2,bp,...
    ux,uy,uz,P,A1,apppsa,atpppsa);
C=Ceval(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
[H,Hiter]=Hcorr(H,X,C,Htol,par,J,dpP0,dpP1,dpP2,bp,...
    ux,uy,uz,P,A1,apppsa,atpppsa);
D=(eye(nq)-X*H*C)*W;
Nu=Nueval(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,...
    apppsa,atpppsa);
qd=D*w+(eye(nq)-X*H*C)*qd0+X*H*Nu;
Gam=Gameval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,...
    uz,P,A1,apppsa,atpppsa);
M=Meval(q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
RM=D'*M*D;
QA=QAEval(t,q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,...
    P,A1,apppsa,atpppsa);
S=Seval(q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
RHS=D'*(M*X*H*Gam+S+QA);
wd=RM\RHS;
Nu=Nueval(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,...
    apppsa,atpppsa);
vd=V'*D*w+V'*(eye(nq)-X*H*C)*qd0+V'*X*H*Nu;
f=[vd;wd];


end



