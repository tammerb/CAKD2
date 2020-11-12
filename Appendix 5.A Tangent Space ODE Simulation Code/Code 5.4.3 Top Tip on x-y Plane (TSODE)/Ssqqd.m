function [Ssq,Ssqd]=Ssqqd(q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);
[rd,pd]=qdPart(qd);

% Evaluate Ssq(q,qd,par) and Ssqd(q,qd,par)
Inertias=[Ixy;Ixy;Iz];
Jpr=diag(Inertias);
Gd=Geval(pd);
G=Geval(p);
T=Teval(Jpr*Gd*p);
Ssq=blkdiag(zeros(3,3),8*Gd'*Jpr*Gd);
Ssqd=blkdiag(zeros(3,3),-8*Gd'*Jpr*G+8*T);

end

