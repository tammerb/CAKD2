function S=SEval(q,qd,par)
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
pd1=[qd(1);qd(2);qd(3);qd(4)];
pd2=[qd(5);qd(6);qd(7);qd(8)];
Gd1=GEval(pd1);
Gd2=GEval(pd2);
S=(16/5)*par(3)*[Gd1'*Gd1*p1;Gd2'*Gd2*p2;0;0;0];


end

