function [qdd,Lam,DE,B]=DAEODEfunct(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

[Pst,Pstt,Pstq,Psttq]=P5(t,q,par);
M=Meval(q,par);
QA=QAeval(t,q,qd,par);
S=Seval(q,qd,par);
Gamma=Gameval(t,q,qd,par);
RHS=[S+QA;-Gamma];
Phisq=P1(t,q,par);
DE=[M,Phisq';Phisq,zeros(nh,nh)];
z=DE\RHS;
i=0;
while i<=nq
i=i+1;
qdd(i)=z(i);
end
j=0;
while j<=nh
Lam(j)=z(nq+j);
end




end

