function P3=P3(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);

% Enter Third Derivative of Constraints P3=((((Psq)qd)sq)qd)sq
[p1,p2,r2] = qpart(q);
[p1d,p2d,r2d] = qdpart(qd);
uz=[0;0;1];
r21=r2+Aeval(p2)*uz+Aeval(p1)*uz;
BT1=BTran(p1,uz);
BT2=BTran(p2,uz);
BT1d=BTran(p1d,uz);
BT2d=BTran(p2d,uz);
e=BT1*p1d+BT2*p2d+r2d;
c=p1d'*BT1'+p2d'*BT2'+r2d';
d=p1d'*BT1d'+p2d'*BT2d';
b1=e'*BT1d+c*BT1d+d*BT1;
b2=e'*BT2d+c*BT2d+d*BT2;
b3=d;
P3=[zeros(2,nq);b1,b2,b3];


end

