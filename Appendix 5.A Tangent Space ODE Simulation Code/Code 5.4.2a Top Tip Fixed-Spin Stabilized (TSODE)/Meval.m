function M=Meval(q,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

% Enter Mass Matrix M=M(q,par)
[r,p]=qPart(q);
G=Geval(p);
inertias=[Ixy,Ixy,Iz];
Jpr=diag(inertias);
M=blkdiag(m*eye(3),4*G'*Jpr*G);

end

