function M=Meval(q,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

I3=eye(3);
p=[q(4);q(5);q(6);q(7)];
Gbar=Gbareval(p);
M=blkdiag(m*I3,4*Gbar'*Jpr*Gbar,zeros(2,2));


end

