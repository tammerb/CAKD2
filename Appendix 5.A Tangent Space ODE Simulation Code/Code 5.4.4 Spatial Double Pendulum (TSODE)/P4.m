function P4=P4(t,q,n,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixx]=parPart(par);
[p1,p2,r2] = qpart(q);
n1=n(1);
n2=n(2);
n3=n(3);
uz=[0;0;1];
r21=r2+Aeval(p1)*uz+Aeval(p2)*uz;
BT1=BTran(p1,uz);
BT2=BTran(p2,uz);

K=Keval(uz,r21);
c1=n1*eye(4)+n3*K+n3*BT1'*BT1;
c2=n2*eye(4)+n3*K+n3*BT2'*BT2;

P4=[c1,n3*BT1'*BT2,n3*BT1';n3*BT2'*BT1,c2,n3*BT2';...
    n3*BT1,n3*BT2,n3*eye(3)];


end

