function [Pst,Pstt,Pstq,Psttq]=P5Eval(tn,q,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

% Enter Constraint t derivatives of P(t,q,qd,par); Default is Zeros

Pst=zeros(nc,1);
Pstt=zeros(nc,1);
Pstq=zeros(nc,ngc);
Psttq=zeros(nc,ngc);

end