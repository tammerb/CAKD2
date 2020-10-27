function [QAsq,QAsqd]=QAsqqdeval(q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);

QAsq=zeros(9,9);
QAsqd=zeros(9,9);


end

