function G=GT(p)
%Evaluate G(p) of Eq.(2.6.2), given p
e0=p(1);
e=[p(2);p(3);p(4)];

G=[-e,-atil(e)+e0*eye(3)];


end
