function P1 = P1(t,q,par)
%Evaluate Constraint Jacobian

[p1,p2,r2] = qpart(q);
uz=[0;0;1];
r21=r2+Aeval(p2)*uz+Aeval(p1)*uz;

BT1=BTran(p1,uz);
BT2=BTran(p2,uz);

P1=[p1',[0,0,0,0],[0,0,0];[0,0,0,0],p2',[0,0,0];r21'*BT1,r21'*BT2,r21'];
end
