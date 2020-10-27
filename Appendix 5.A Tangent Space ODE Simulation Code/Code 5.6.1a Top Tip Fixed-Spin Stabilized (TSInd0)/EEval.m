function Ebar=EEval(p)
e0=p(1);
e=[p(2);p(3);p(4)];
Ebar=[-e,atil(e)+e0*eye(3)];


end

