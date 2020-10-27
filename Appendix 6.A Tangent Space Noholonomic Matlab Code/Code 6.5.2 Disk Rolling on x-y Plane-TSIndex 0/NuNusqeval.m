function [Nu,Nusq]=NuNusqeval(t,q,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);

Nu=zeros(nh+nd,1);
Nusq=zeros(nh+nd,nq);


end

