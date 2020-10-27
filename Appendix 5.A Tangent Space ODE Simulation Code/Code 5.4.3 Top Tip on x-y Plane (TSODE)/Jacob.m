function [Rv,Rvd,Rvdd]=Jacob(t,q,qd,qdd,U,V,B,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

Phisq=P1(t,q,par);
QA=QAeval(t,q,qd,par);
[QAsq,QAsqd]=QAsqqd(t,q,qd,par);
S=Seval(q,qd,par);
[Ssq,Ssqd]=Ssqqd(q,qd,par);
Gamma=Gameval(t,q,qd,par);
[Gammasq,Gammasqd]=Gammasqqd(t,q,qd,par);
M=Meval(q,par);
vd=V'*qd;
vdd=V'*qdd;
[B,Biter]=Beval(t,q,B,U,par);
D=Deval(t,q,U,V,B,par);

Rvdd=D'*M*D;
Rvd=-D'*(M*U*B*Gammasqd+Ssqd+QAsqd)*D;
Rv=D'*(-M*U*B*P2(t,q,2*D*vdd-U*B*Gamma,par)+...
    P4(t,q,B'*U'*(M*U*B*Gamma+S+QA),par)+...
    M2(q,D*vdd-U*B*Gamma,par)-M*U*B*Gammasq-Ssq-QAsq-...
    (M*U*B*Gammasqd+Ssqd+QAsqd)*U*B*(P2(t,q,-D*vd,par)))*D;



end

