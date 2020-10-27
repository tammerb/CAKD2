function C=Ceval(t,q,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

C=[Phiqeval(t,q,par);Eeval(t,q,par)];

end

