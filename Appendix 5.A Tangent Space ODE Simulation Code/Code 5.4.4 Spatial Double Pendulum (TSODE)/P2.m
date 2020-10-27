function P2 = P2(t,q,x,par)

[p1,p2,r2] = qpart(q);
[xp1,xp2,xr2] = qdpart(x);
uz=[0;0;1];
r21=r2+Aeval(p2)*uz+Aeval(p1)*uz;

BT1=BTran(p1,uz);
BT2=BTran(p2,uz);

a=BT1*xp1+BT2*xp2+xr2;

P2=[xp1',[0,0,0,0],[0,0,0];[0,0,0,0],xp2',[0,0,0];...
    a'*BT1+r21'*BTran(xp1,uz),a'*BT2+r21'*BTran(xp2,uz),a'];

end

