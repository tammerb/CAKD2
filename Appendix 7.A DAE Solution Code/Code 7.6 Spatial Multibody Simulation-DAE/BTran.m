function BT=BTran(p,apr)
e0=p(1);
e=[p(2);p(3);p(4)];
I3=eye(3);
etil=atil(e);
BT=2*[(e0*I3+etil)*apr,e*apr'-(e0*I3+etil)*atil(apr)];


end

