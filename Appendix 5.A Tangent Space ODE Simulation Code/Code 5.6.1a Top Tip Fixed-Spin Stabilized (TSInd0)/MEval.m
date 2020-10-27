function M=MEval(q,par)
Inertias=[par(3);par(3);par(4)];
Jpr=diag(Inertias);
pq=[q(4);q(5);q(6);q(7)];
G=GEval(pq);
M=blkdiag(par(1)*eye(3),4*G'*Jpr*G);


end

