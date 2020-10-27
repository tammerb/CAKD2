function M=Meval(q,par,Jpr)
%Enter Mass Matrix nqxnq; example below is for rolling eooipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
I3=eye(3);
p=[q(4);q(5);q(6);q(7)];
G=Geval(p);
M=blkdiag(m*I3,4*G'*Jpr*G,zeros(5,5));


end

