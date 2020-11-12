function S=SEval(q,qd,par)
Inertias=[par(3);par(4);par(5)];
Jpr=diag(Inertias);
pq=[q(4);q(5);q(6);q(7)];
pqd=[qd(4);qd(5);qd(6);qd(7)];
Gd=GEval(pqd);
S=[0;0;0;8*Gd'*Jpr*Gd*pq];


end

