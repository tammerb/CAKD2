function Gam=GamEval(tn,q,qd,SJDT,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

[P,Pst,Pstt]=P5Eval(tn,q,SJDT,par);
P2=P2Eval(q,qd,SJDT,par);
Gam=P2*qd+Pstt;

end

