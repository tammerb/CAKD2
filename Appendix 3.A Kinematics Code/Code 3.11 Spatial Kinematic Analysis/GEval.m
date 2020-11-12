function G=GEval(x)
e0=x(1);
e=[x(2);x(3);x(4)];
etil=atil(e);
G=[-e,-etil+e0*eye(3)];

end

