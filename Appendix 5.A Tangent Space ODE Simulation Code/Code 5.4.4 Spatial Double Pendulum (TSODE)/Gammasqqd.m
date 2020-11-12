function [Gammasq,Gammasqd]=Gammasqqd(t,q,qd,par)

[Pst,Pstt,Pstq,Psttq]=P5(t,q,par);
Pstqmqdsq=P6(t,q,qd,par);
Gammasq=P3(t,q,qd,par)+2*Pstqmqdsq+Psttq;
Gammasqd=2*P2(t,q,qd,par)+2*Pstq;


end

