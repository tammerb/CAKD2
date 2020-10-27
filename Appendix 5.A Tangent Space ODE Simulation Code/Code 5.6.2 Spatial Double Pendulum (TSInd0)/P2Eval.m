function P2=P2Eval(q,x)
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
r2=[q(9);q(10);q(11)];
x1=[x(1);x(2);x(3);x(4)];
x2=[x(5);x(6);x(7);x(8)];
x3=[x(9);x(10);x(11)];
kpr=[0;0;1];
r21=r2+ATran(p2)*kpr+ATran(p1)*kpr;
a=BTran(p1,kpr)*x1+BTran(p2,kpr)*x2+x3;
P2=[x1',[0,0,0,0],[0,0,0];[0,0,0,0],x2',[0,0,0];...
    a'*BTran(p1,kpr)+r21'*BTran(x1,kpr),...
    a'*BTran(p2,kpr)+r21'*BTran(x2,kpr),a'];
end

