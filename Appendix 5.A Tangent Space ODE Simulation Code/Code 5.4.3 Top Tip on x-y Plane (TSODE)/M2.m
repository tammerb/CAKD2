function M2=M2(q,mu,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

% Enter M2=(M(q,par)a)sq
p=[q(4);q(5);q(6);q(7)];
mup=[mu(4);mu(5);mu(6);mu(7)];
G=Geval(p);
Jpr=diag([Ixy;Ixy;Iz]);
M2=[zeros(3,7);zeros(4,3),4*Teval(Jpr*G*mup)-4*G'*Jpr*Geval(mup)];


end

