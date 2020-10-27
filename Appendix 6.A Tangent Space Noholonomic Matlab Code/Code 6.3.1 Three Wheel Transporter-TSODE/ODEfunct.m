function [f,RM]=ODEfunct(t,y,ue,V,W,X,q0,qd0,par)
%Example is transporter
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);

[v,w] = yPart(y)

q=q0+V*v;
phi=q(3);

C=Ceval(t,q,par);
a=C*X;
D=(eye(nq)-(1/a)*X*C)*W;
qd=D*w+(eye(nq)-(1/a)*X*C)*qd0;
phid=qd(3);
Gam=Gameval(t,q,qd,par);
M=Meval(q,par);
RM=D'*M*D;
QA=QAEval(t,q,qd,par);
AT=ATran(phi);
spc=[0;1];
S=-[m*(phid^2)*AT*spc;0];
RHS=D'*QA-(1/a)*D'*M*X*Gam-D'*S;
wd=RM\RHS;
vd=D*w+(eye(nq)-(1/a)*X*C)*qd0;
f=[vd;wd];


end



