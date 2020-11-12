function C2=C2eval(t,q,x,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

P2=P2eval(t,q,x,par);
E2=E2eval(q,x,par);

C2=[P2;E2];


end

