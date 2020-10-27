function [Gam,Gamsq,Gamsqd] = GamEval(q,qd,par)

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

P2=P2Eval(q,qd);
Gam=P2*qd;
Gamsq=[zeros(nh,nq)];
Gamsqd=2*P2;

end

