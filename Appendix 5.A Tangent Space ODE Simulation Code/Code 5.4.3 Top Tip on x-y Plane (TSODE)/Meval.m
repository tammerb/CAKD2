function M=Meval(q,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

% Enter Mass Matrix M=M(q,par)
[r,p]=qPart(q);
Inertias=[Ixy;Ixy;Iz];
Jpr=diag(Inertias);
G=Geval(p);
M=blkdiag(m*eye(3),4*G'*Jpr*G);

end

