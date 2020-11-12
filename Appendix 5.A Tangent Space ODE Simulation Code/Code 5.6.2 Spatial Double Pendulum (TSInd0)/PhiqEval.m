function Phiq = PhiqEval(q)
%Evaluate Constraint Jacobian
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
r2=[q(9);q(10);q(11)];
AT1=ATran(p1);
AT2=ATran(p2);
kpr=[0;0;1];
r21=r2+AT1*kpr+AT2*kpr;
BT1=BTran(p1,kpr);
BT2=BTran(p2,kpr);
Phiq=[p1',0,0,0,0,0,0,0;0,0,0,0,p2',0,0,0;r21'*BT1,r21'*BT2,r21'];
    end
