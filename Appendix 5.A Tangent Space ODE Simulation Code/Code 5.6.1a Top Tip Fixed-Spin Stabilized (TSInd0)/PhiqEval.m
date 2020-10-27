function Phiq = PhiqEval(q)
%Evaluate Constraint Jacobian
p=[q(4);q(5);q(6);q(7)];
I3=eye(3);
kpr=[0;0;1];
BT=BTran(p,kpr);
Phiq=[I3,-BT;[0,0,0],p'];
    end
