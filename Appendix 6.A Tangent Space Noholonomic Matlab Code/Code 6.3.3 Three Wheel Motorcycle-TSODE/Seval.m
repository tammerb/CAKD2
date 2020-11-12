function S=Seval(q,qd,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa)
%Quadratic velocity term is nqx1; example below is for rolling ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
[r,p,a,s]=qPart(q);
[rd,pd,ad,sd]=qdPart(qd);
Gd=Geval(pd);
S=[0;0;0;8*Gd'*J*Gd*p;0;0;0];


end

