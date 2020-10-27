function [H,Hiter]=Hcorr(H,X,C,Htol,par,J,dpP0,dpP1,dpP2,bp,...
    ux,uy,uz,P,A1,apppsa,atpppsa)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
i=1;
Herrnorm=Htol+1;
I=eye(nh+nd);
while Herrnorm>Htol
i=i+1;
delH=-H*C*X*H+H;
H=H+delH;
Herr=C*X*H-I;
Herrnorm=norm(Herr);
Hiter=i-1;
end

