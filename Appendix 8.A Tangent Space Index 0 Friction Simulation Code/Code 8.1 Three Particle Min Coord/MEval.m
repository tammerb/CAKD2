function M=MEval(q,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

m=[m1;m2;m3];

M=diag(m);


end

