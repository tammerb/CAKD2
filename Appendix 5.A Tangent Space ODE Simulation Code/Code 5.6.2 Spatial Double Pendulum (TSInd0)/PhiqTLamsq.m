function PqTLsq = PhiqTLamsq(q,Lam)
p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
r2=[q(9);q(10);q(11)];
AT1=ATran(p1);
AT2=ATran(p2);
kpr=[0;0;1];
r21=r2+AT1*kpr+AT2*kpr;
kpr=[0;0;1];
BT1=BTran(p1,kpr);
BT2=BTran(p2,kpr);
K=KEval(kpr,r21);
c1=Lam(1)*eye(4)+Lam(3)*K+Lam(3)*BT1'*BT1;
c2=Lam(2)*eye(4)+Lam(3)*K+Lam(3)*BT2'*BT2;
PqTLsq=[c1,Lam(3)*BT1'*BT2,Lam(3)*BT1';Lam(3)*BT2'*BT1,c2,Lam(3)*BT2';...
    Lam(3)*BT1,Lam(3)*BT2,Lam(3)*eye(3)];


end

