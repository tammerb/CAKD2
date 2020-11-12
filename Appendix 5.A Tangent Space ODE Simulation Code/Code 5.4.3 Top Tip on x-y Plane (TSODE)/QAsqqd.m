function [QAsq,QAsqd]=QAsqqd(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);
[rd,pd]=qdPart(qd);

% Evaluate Derivatives QAsq(t,q,qd,par) and QAsqd(t,q,qd,par)

AT=ATran(p);
uz=[0;0;1];
BT=BTran(p,uz);
BTd=BTran(pd,uz);
Kev = Keval(uz,K*(r-AT*uz)+C*(rd-BT*pd));
QAsq=[-K*eye(3),K*BT+C*BTd;K*BT',Kev-BT'*(K*BT+C*BTd)];
QAsqd=C*[-eye(3),BT;BT',-BT'*BT];   

end

