function M2=M2(q,mu,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);
[p1,p2,r2] = qpart(q);
[mup1,mup2,mur2] = qpart(mu);

% Enter M2=(M(q,par)a)sq
Gp1=Geval(p1);
Gp2=Geval(p2);
Gmup1=Geval(mup1);
Gmup2=Geval(mup2);
T1=Teval(Gp1*mup1);
T2=Teval(Gp2*mup2);
M2=(8*m/5)*blkdiag(-Gp1'*Gmup1+T1,-Gp2'*Gmup2+T2,zeros(3,3));


end

