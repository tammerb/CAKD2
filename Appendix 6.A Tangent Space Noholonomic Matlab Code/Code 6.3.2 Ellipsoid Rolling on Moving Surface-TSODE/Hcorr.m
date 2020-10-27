function [H,Hiter]= Hcorr(H,X,C,Htol,par)
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
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

