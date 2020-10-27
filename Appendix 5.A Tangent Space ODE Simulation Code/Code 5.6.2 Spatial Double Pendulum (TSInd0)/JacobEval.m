function [Rv,Rvd,Rvdd]=JacobEval(q,qd,vdd,Lam,U,V,B,par)

[Sq,Sqd]=SqqdEval(q,qd,par);
M=MEval(q,par);
PqTLsq=PhiqTLamsq(q,Lam);
[Gam,Gamsq,Gamsqd] = GamEval(q,qd);
qdd=V*vdd-U*B*Gam;

Rvdd=M*V;
Rvd=-(M*U*B*Gamsqd+Sqd)*V;
Rv=(M2Eval(q,qdd,par)+M*(-U*B*(P2Eval(q,qdd)+Gamsq)...
    +U*B*Gamsqd*U*B*P2Eval(q,qd))+PqTLsq-Sq+Sqd*U*B*P2Eval(q,qd))*V;


end



