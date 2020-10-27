function S=Seval(q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);
[p1,p2,r2] = qpart(q);
[pd1,pd2,rd2] = qdpart(qd);

% Evaluate S=S(q,qd,par)
Gd1=Geval(pd1);
Gd2=Geval(pd2);
S=(16/5)*m*[Gd1'*Gd1*p1;Gd2'*Gd2*p2;0;0;0];

end

