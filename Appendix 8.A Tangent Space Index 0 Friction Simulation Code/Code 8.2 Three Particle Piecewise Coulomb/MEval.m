function M=MEval(q,par)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

m=[m1;m2;m3];

M=diag(m);


end

