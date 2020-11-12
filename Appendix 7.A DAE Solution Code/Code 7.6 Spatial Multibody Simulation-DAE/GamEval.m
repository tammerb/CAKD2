function Gam=GamEval(tn,q,qd,SJDT,par)


[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
P2=P2Eval(tn,q,qd,SJDT,par);
Gam=P2*qd+Pstt;

end

