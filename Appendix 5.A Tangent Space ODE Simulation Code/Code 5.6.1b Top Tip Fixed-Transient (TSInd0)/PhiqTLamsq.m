function PqTLsq = PhiqTLamsq(Lam)
kpr=[0;0;1];
Lamk=[Lam(1);Lam(2);Lam(3)];
Lamp=Lam(4);
PqTLsq=[zeros(3,7);zeros(4,3),-KEval(kpr,Lamk)+Lamp*eye(4)];


end

