function [Gam,Gamsq,Gamsqd] = GamEval(q,qd)
P2=P2Eval(q,qd);
P3=P3Eval(q,qd);
Gam=P2*qd;
Gamsq=P3;
Gamsqd=2*P2;

end

