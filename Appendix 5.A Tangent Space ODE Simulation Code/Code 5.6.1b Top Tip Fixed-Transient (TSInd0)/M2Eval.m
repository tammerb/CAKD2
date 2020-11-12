function M2=M2Eval(q,x,par)

Inertias=[par(3);par(4);par(5)];
Jpr=diag(Inertias);
pq=[q(4);q(5);q(6);q(7)];
px=[x(4);x(5);x(6);x(7)];
G=GEval(pq);
a=Jpr*G*px;
T=TEval(a);
M2=blkdiag(zeros(3,3),-4*G'*Jpr*GEval(px)+4*T);


end

