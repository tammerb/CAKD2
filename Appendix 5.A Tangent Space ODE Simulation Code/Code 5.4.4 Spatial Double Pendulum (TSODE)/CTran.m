function CT=CTran(p,a)
e0=p(1);
e=[p(2);p(3);p(4)];
I3=eye(3);
CT=2*[(e0*I3-atileval(e))*a,e*a'+(e0*I3-atileval(e))*atileval(a)];


end

