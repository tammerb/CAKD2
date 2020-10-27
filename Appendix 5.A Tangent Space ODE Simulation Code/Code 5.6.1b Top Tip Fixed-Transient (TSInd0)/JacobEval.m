function [Rv,Rvd,Rvdd]=JacobEval(q,qd,vdd,Lam,U,V,B,par)

Phiq = PhiqEval(q);
P20=P2Eval(q,qd);
QA=QAEval(q,qd,par);
S=SEval(q,qd,par);
[Sq,Sqd]=SqqdEval(q,qd,par);
M=MEval(q,par);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
PqTLsq = PhiqTLamsq(Lam);

Rvdd=M*V;
Rvd=-(M*U*B*Gamsqd+Sqd)*V;
Rv=(M2Eval(q,V*vdd,par)-M*U*B*P2Eval(q,V*vdd)-M2Eval(q,U*B*Gam,par)...
    +M*U*B*P2Eval(q,U*B*Gam)-M*U*B*Gamsq+PqTLsq-Sq...
    +(M*U*B*Gamsqd+Sqd)*U*B*P20)*V;


end



