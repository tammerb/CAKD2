function Gam=GamEval(tn,q,qd,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
P2=P2Eval(tn,q,qd,SJDT,par);
Gam=P2*qd+Pstt;

end

