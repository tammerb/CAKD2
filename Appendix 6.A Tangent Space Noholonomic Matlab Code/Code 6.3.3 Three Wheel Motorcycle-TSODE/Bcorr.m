function [B,Biter]=Bcorr(t,q,B,U,Btol,par,J,dpP0,dpP1,dpP2,bp,...
    ux,uy,uz,P,A1,apppsa,atpppsa)
i=1;
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
Berrnorm=Btol+1;
I=eye(nu);
Phiq=P1(t,q,par,J,dpP0,dpP1,dpP2,bp,ux,uy,uz,P,A1,apppsa,atpppsa);
while Berrnorm>Btol
i=i+1;
delB=-B*Phiq*U*B+B;
B=B+delB;
Berr=B*Phiq*U-I;
Berrnorm=norm(Berr);
Biter=i-1;
end

end

