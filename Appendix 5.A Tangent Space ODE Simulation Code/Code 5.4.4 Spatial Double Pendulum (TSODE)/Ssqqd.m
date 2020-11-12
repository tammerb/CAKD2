function [Ssq,Ssqd]=Ssqqd(q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);

% Evaluate Ssq(q,qd,par) and Ssqd(q,qd,par)
[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);
[p1,p2,r2] = qpart(q);
[p1d,p2d,r2d] = qdpart(qd);

G1=Geval(p1);
G2=Geval(p2);
G1d=Geval(p1d);
G2d=Geval(p2d);
Ssq=(16/5)*m*blkdiag(G1d'*G1d,G2d'*G2d,zeros(3,3));
Ssqd=(16/5)*m*blkdiag(-G1d'*G1+Teval(G1d*p1),-G2d'*G2+Teval(G2d*p2),...
    zeros(3,3));

end

