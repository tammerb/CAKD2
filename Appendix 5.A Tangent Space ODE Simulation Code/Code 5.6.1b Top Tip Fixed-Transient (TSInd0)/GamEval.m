function [Gam,Gamsq,Gamsqd] = GamEval(q,qd)
P2=P2Eval(q,qd);
Gam=P2*qd;
Gamsq=[zeros(4,7)];
Gamsqd=2*P2;

end

