function G=Geval(x)
e0=x(1);
e=[x(2);x(3);x(4)];
etil=[0,-x(4),x(3);x(4),0,-x(2);-x(3),x(2),0];
G=[-e,-etil+e0*eye(3)];


end

