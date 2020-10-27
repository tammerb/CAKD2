function [B,Biter]=Bcorr(t,q,B,U,par)
i=1;
[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
Berrnorm=Btol+1;
I=eye(nu);
Phiq=Phiqeval(t,q,par);
while Berrnorm>Btol
i=i+1;
delB=-B*Phiq*U*B+B;
B=B+delB;
Berr=B*Phiq*U-I;
Berrnorm=norm(Berr);
Biter=i-1;
end

end

