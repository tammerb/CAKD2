function A=BAeval(p)
e0=p(1);
e=[p(2);p(3);p(4)];
I3=eye(3);
etil=Batil(e);
A=(e0^2-e'*e)*I3+2*e*e'+2*e0*etil;


end
