function S=Seval(q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

% Evaluate S=S(q,qd,par)
[r,p]=qPart(q);
[rd,pd]=qdPart(qd);
Gd=Geval(pd);
inertias=[Ixy,Ixy,Iz];
Jpr=diag(inertias);

S=[zeros(3,1);8*Gd'*Jpr*Gd*p];

end

