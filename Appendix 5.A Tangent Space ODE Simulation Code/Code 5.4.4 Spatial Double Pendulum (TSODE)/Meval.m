function M=Meval(q,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);
[p1,p2,r2] = qpart(q);

% Enter Mass Matrix M=M(q,par)
G1=Geval(p1);
G2=Geval(p2);
M=m*blkdiag((8/5)*G1'*G1,(8/5)*G2'*G2,eye(3));

end

