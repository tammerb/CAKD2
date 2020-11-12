function [Sq,Sqd]=SqqdEval(q,qd,par)
c=[par(3);par(3);par(4)];
Jpr=diag(c);
pq=[q(4);q(5);q(6);q(7)];
pqd=[qd(4);qd(5);qd(6);qd(7)];
Gd=GEval(pqd);
G=GEval(pq);
a1=Jpr*Gd*pq;
Sq=blkdiag(0*eye(3),8*Gd'*Jpr*Gd);
Sqd=blkdiag(0*eye(3),-8*Gd'*Jpr*G+8*TEval(a1));


end

