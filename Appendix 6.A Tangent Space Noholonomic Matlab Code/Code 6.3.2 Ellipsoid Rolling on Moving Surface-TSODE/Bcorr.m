function [B,Biter]=Bcorr(t,q,B,X,Btol,par,L)
i=1;
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
Berrnorm=Btol+1;
I=eye(nu);
Phiq=Phiqeval(t,q,par,L);
while Berrnorm>Btol
i=i+1;
delB=-B*Phiq*X*B+B;
B=B+delB;
Berr=B*Phiq*X-I;
Berrnorm=norm(Berr);
Biter=i-1;
end

end

