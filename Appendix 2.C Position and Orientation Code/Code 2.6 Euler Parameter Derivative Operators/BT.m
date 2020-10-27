function BT=BT(p,apr)
%Evaluate B(p,apr) of Eq.(2.6.25), given p and apr
e0=p(1);
e=[p(2);p(3);p(4)];
I3=eye(3);
etil=atil(e);
BT=2*[(e0*I3+etil)*apr,e*apr'-(e0*I3+etil)*atil(apr)];

end

