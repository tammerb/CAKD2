function [B,Biter]=Bcorr(t,q,B,X,Btol,par)
i=1;
[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ]=Partpar(par);
Berrnorm=Btol+1;
I=eye(nu);
Phiq = P1(t,q,par);
while Berrnorm>Btol
i=i+1;
delB=-B*Phiq*X*B+B;
B=B+delB;
Berr=B*Phiq*X-I;
Berrnorm=norm(Berr);
Biter=i-1;
end

end

