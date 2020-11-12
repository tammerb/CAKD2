function P3=P3Eval(q,qd)
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
r2=[q(9);q(10);q(11)];
p1d=[qd(1);qd(2);qd(3);qd(4)];
p2d=[qd(5);qd(6);qd(7);qd(8)];
r2d=[qd(9);qd(10);qd(11)];
kpr=[0;0;1];
a=BTran(p1,kpr)*p1d+BTran(p2,kpr)*p2d+r2d;
c=p1d'*BTran(p1,kpr)'+p2d'*BTran(p2,kpr)'+r2d';
d=p1d'*BTran(p1d,kpr)'+p2d'*BTran(p2d,kpr)';
b1=(a'+c)*BTran(p1d,kpr)+d*BTran(p1,kpr);
b2=(a'+c)*BTran(p2d,kpr)+d*BTran(p2,kpr);
b3=d;
P3=[zeros(2,11);b1,b2,b3];


end

