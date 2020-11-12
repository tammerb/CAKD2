function S=Seval(q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);
[rd,pd]=qdPart(qd);

% Evaluate S=S(q,qd,par)
Inertias=[Ixy;Ixy;Iz];
Jpr=diag(Inertias);

Gd=Geval(pd);
S=[0;0;0;8*Gd'*Jpr*Gd*p];

end

