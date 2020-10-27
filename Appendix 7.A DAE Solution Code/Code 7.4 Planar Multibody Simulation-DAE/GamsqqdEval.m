function [Gamsq,Gamsqd]=GamsqqdEval(tn,q,qd,PJDT,par)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA,NRSDA]=parPart(par);

P3=P3Eval(tn,q,qd,PJDT,par);
P2=P2Eval(tn,q,qd,PJDT,par);
Gam=P2*qd;
Gamsq=P3;
Gamsqd=2*P2;

end

