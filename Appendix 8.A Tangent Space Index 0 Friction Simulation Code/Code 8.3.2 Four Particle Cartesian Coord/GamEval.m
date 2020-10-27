function [Gam,Gamsq,Gamsqd] = GamEval(q,qd,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

P2=P2Eval(q,qd);
Gam=P2*qd;
Gamsq=[zeros(nh,nq)];
Gamsqd=2*P2;

end

