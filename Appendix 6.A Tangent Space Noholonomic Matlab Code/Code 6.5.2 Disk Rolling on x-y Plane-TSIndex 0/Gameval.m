function Gam=Gameval(t,q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[rd,pd,ad]=qdPart(qd);

Gam=[P2eval(t,q,qd,par)*qd;E2eval(q,qd,par)*qd];

end



