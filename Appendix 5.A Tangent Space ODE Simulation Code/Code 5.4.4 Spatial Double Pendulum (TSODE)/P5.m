function [Pst,Pstt,Pstq,Psttq]=P5(t,q,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);

% Enter Constraint t derivatives of P(t,q,qd,par); Default is Zeros

Pst=zeros(nh,1);
Pstt=zeros(nh,1);
Pstq=zeros(nh,nq);
Psttq=zeros(nh,nq);


end

