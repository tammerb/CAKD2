function M2=M2(q,mu,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

% Enter M2=(M(q,par)a)sq
[r,p]=qPart(q);
[mur,mup]=qPart(mu);
G=Geval(p);
Gmup=Geval(mup);
inertias=[Ixy,Ixy,Iz];
Jpr=diag(inertias);

M2=[zeros(3,nq);zeros(4,3),-4*G'*Jpr*Gmup+4*Teval(Jpr*G*mup)];


end

