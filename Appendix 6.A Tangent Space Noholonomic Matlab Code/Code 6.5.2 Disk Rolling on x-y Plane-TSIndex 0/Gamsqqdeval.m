function [Gamsq,Gamsqd] = Gamsqqdeval(t,q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

Gamsq=[P3eval(t,q,qd,par);E3eval(q,qd,par)];
Gamsqd=[2*P2eval(t,q,qd,par);2*E2eval(q,qd,par)];


end

