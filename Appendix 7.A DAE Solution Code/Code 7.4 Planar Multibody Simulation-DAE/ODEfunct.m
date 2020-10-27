function [qdd,Lam,ECond]=ODEfunct(tn,q,qd,PMDT,PTSDAT,PRSDAT,PJDT,par)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA,NRSDA]=parPart(par);

Phiq=PhiqEval(tn,q,PJDT,par);
M=MEval(PMDT,par);
E=[M,Phiq';Phiq,zeros(nc,nc)];
Gam=GamEval(tn,q,qd,PJDT,par);
QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par);
RHS=[QA;-Gam];

x=E\RHS;

Pqdd=[eye(ngc),zeros(ngc,nc)];
PLam=[zeros(nc,ngc),eye(nc)];
qdd=Pqdd*x;
Lam=PLam*x;
ECond=cond(E);

end




