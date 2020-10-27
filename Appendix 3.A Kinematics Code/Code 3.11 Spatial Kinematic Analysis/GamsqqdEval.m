function [Gamsq,Gamsqd]=GamsqqdEval(tn,q,qd,SJDT,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

P2=P2Eval(tn,q,qd,SJDT,par);
P3=P3Eval(tn,q,qd,SJDT,par);
Gamsq=P3;
Gamsqd=2*P2;

