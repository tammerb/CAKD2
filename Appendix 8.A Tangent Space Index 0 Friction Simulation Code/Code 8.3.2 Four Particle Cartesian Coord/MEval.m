function M=MEval(q,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

m=[m1;m1;m1;m2;m2;m2;m3;m3;m3;m4];

M=diag(m);


end

