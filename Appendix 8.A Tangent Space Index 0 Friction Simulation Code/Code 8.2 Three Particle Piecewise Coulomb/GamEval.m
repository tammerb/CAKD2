function [Gam,Gamsq,Gamsqd] = GamEval(q,qd,par,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);


P2=P2Eval(q,qd,par,mode);
Gam=P2*qd;

if mode==1
Gamsq=[zeros(nh1,nq)];
end
if mode>1
Gamsq=[zeros(nh2,nq)];
end

Gamsqd=2*P2;


end

