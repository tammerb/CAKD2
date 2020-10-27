function [B,Biter]=Bcorr(t,q,B,U,Btol,par)
i=1;
[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ,h,utol,Btol,Htol]=Partpar(par);
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

