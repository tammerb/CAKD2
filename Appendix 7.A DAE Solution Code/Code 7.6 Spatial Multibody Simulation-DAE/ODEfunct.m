function [qdd,Lam,ECond]=ODEfunct(tn,q,qd,SMDT,STSDAT,SJDT,par)

[nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA]=parPart(par);

Gam=GamEval(tn,q,qd,SJDT,par);
QA=QAEval(tn,q,qd,SMDT,STSDAT,par);
S=SEval(q,qd,SMDT,par);
RHS=[QA+S;-Gam];
M=MEval(q,SMDT,par);
Phiq=PhiqEval(tn,q,SJDT,par);
E=[M,Phiq';Phiq,zeros(nc,nc)];
ECond=cond(E);

x=E\RHS;

qdd=zeros(ngc,1);
i=1;
while i<=ngc
qddi=x(i);
qdd=Add(qdd,qddi,i-1,0);
i=i+1;
end
Lam=zeros(nc,1);
i=1;
while i<=nc
Lami=x(ngc+i);
Lam=Add(Lam,Lami,i-1,0);
i=i+1;
end

end




