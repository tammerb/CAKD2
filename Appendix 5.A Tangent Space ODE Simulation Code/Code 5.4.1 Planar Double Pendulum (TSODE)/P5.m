function [Pst,Pstt,Pstq,Psttq]=P5(t,q,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

% Enter Constraint t derivatives of P(t,q,qd,par); Default is Zeros

Pst=zeros(nh,1);
Pstt=zeros(nh,1);
Pstq=zeros(nh,nq);
Psttq=zeros(nh,nq);


end

