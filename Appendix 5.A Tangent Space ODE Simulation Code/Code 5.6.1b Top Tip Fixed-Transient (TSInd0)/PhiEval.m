function Phi = PhiEval(q)
%Evaluate constraint function
r=[q(1);q(2);q(3)];
p=[q(4);q(5);q(6);q(7)];
kpr=[0;0;1];
AT=ATran(p);
Phi=[r-AT*kpr;(p'*p-1)/2];
end

