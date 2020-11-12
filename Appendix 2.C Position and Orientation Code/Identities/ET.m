function E=ET(p)
%Evaluate E(p) of Eq.(2.6.1), given p 
e0=p(1);
e=[p(2);p(3);p(4)];

E=[-e,atil(e)+e0*eye(3)];

end

