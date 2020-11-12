function G=Geval(p)
e0=p(1);
e=[p(2);p(3);p(4)];
etil=[0,-p(4),p(3);p(4),0,-p(2);-p(3),p(2),0];
G=[-e,-etil+e0*eye(3)];


