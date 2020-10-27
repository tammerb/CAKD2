function [q0,qd0,V,U,W,X,B,H,jRepar]=Param(t,q,qd,jRepar,par,J,...
    dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
jRepar=jRepar+1;
q0=q;
qd0=qd;
Phiq=P1(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
V=null(Phiq);
U=Phiq';
B=inv(U'*U);
C0=Ceval(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
X=C0';
W=null(C0);
H=inv(C0*C0');

end

