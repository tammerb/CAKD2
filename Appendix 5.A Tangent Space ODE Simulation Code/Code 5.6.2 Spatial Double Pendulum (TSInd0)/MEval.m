function M=MEval(q,par)
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
G1=GEval(p1);
G2=GEval(p2);
M=par(3)*blkdiag(4*0.4*G1'*G1,4*0.4*G2'*G2,eye(3));


end

