function Gam=GamEval(tn,q,qd,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

P2=P2Eval(tn,q,qd,PJDT,par);
[Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par);
Gam=P2*qd+Pstt;

end

