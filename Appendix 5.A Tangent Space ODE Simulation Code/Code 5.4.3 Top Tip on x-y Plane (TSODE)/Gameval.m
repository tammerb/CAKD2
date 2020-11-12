function Gamma=Gameval(t,q,qd,par)

[Pst,Pstt,Pstq,Psttq]=P5(t,q,par);
Gamma=P2(t,q,qd,par)*qd+2*Pstq*qd+Pstt;


end

