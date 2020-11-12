function P2 = P2(t,q,x,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);
[r,p]=qPart(q);
xp=[x(4);x(5);x(6);x(7)];

uz=[0;0;1];
BT=BTran(xp,uz);
P2=[0,0,0,-uz'*BT;0,0,0,xp'];
end

