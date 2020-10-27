function M2=M2eval(q,x,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[xr,xp,xa]=xPart(x);


M2=[zeros(3,9);zeros(4,3),4*Tbareval(Jpr*Gbareval(p)*xp)-...
    4*Gbareval(p)'*Jpr*Gbareval(xp),zeros(4,2);zeros(2,9)];


end

