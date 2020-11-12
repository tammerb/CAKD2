function Phi = PhiEval(q)
%Evaluate constraint function
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
r2=[q(9);q(10);q(11)];
kpr=[0;0;1];
AT1=ATran(p1);
AT2=ATran(p2);
r21=r2+AT1*kpr+AT2*kpr;
Phi=[(p1'*p1-1)/2;(p2'*p2-1)/2;(r21'*r21-1)/2];
end

