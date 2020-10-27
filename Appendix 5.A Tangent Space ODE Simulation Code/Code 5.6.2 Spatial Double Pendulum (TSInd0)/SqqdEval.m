function [Sq,Sqd]=SqqdEval(q,qd,par)
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
p1d=[qd(1);qd(2);qd(3);qd(4)];
p2d=[qd(5);qd(6);qd(7);qd(8)];
G1=GEval(p1);
G2=GEval(p2);
G1d=GEval(p1d);
G2d=GEval(p2d);
Sq=(16/5)*par(3)*blkdiag(G1d'*G1d,G2d'*G2d,0*eye(3));
Sqd=(16/5)*par(3)*blkdiag(-G1d'*G1+TEval(G1d*p1),-G2d'*G2+TEval(G2d*p2),...
    zeros(3,3));


end

