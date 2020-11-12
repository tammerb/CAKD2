function [Ssq,Ssqd]=Ssqqd(q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

% Evaluate Ssq(q,qd,par) and Ssqd(q,qd,par)
[r,p]=qPart(q);
[rd,pd]=qdPart(qd);
G=Geval(p);
Gd=Geval(pd);
inertias=[Ixy,Ixy,Iz];
Jpr=diag(inertias);

Ssq=[zeros(3,nq);zeros(4,3),8*Gd'*Jpr*Gd];
Ssqd=[zeros(3,nq);zeros(4,3),-8*Gd'*Jpr*G+8*Teval(Jpr*Gd*p)];

end

