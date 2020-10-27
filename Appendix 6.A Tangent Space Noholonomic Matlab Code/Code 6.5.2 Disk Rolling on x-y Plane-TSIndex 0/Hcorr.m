function [H,Hiter]= Hcorr(H,X,C,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

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

