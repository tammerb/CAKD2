function CT=CT(p,a)
%Evaluate C(p,apr) of Eq.(2.6.31), given p and apr
e0=p(1);
e=[p(2);p(3);p(4)];
I3=eye(3);
CT=2*[(e0*I3-atil(e))*a,e*a'+(e0*I3-atil(e))*atil(a)];


end

