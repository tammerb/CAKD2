function QA=QAeval(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);
[rd,pd]=qdPart(qd);

uz=[0;0;1];
BT=BTran(p,uz);
rt=r-ATran(p)*uz;
rtd=rd-BT*pd;

QA=[-m*g*uz-K*rt-C*rtd;BT'*(K*rt+C*rtd)];
end

