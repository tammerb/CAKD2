function P2=P2Eval(q,x)
xp=[x(4);x(5);x(6);x(7)];
I3=eye(3);
kpr=[0;0;1];
BT=BTran(xp,kpr);
P2=[zeros(3,3),-BT;0,0,0,xp'];
end

