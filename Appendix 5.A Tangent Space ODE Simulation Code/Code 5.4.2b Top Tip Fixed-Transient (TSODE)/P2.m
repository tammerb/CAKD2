function P2 = P2(t,q,x,par)

[r,p]=qPart(q);
[xr,xp]=qPart(x);
uz=[0;0;1];

P2=[zeros(3,3),-BTran(xp,uz);0,0,0,xp'];
end

