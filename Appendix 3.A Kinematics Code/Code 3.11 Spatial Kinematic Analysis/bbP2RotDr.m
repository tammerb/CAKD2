function [P21,P22]=bbP2RotDr(i,j,vx1pr,vy1pr,vx2pr,q,qd,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

[r1,p1]=qPart(q,i);
[xr1,xp1]=qPart(qd,i);
A1=ATran(p1);
BT1y=BTran(p1,vy1pr);
BT1x=BTran(p1,vx1pr);
BTx1y=BTran(xp1,vy1pr);
BTx1x=BTran(xp1,vx1pr);

if j==0
c=vx1pr'*A1'*vx2pr;
s=vy1pr'*A1'*vx2pr;
if abs(c)>=abs(s)
P21=[0,0,0,vx2pr'*BTx1y];
else
P21=[0,0,0,vx2pr'*BTx1x];
end
P22=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
[xr2,xp2]=qPart(qd,j);
A2=ATran(p2);
BT2x=BTran(p2,vx2pr);
BTx2x=BTran(xp2,vx2pr);
c=vx1pr'*A1'*A2*vx2pr;
s=vy1pr'*A1'*A2*vx2pr;
if abs(c)>=abs(s)
P21=[0,0,0,vx2pr'*A2'*BTx1y+xp2'*BT2x'*BT1y];
P22=[0,0,0,vy1pr'*A1'*BTx2x+xp1'*BT1y'*BT2x];
else
P21=[0,0,0,vx2pr'*A2'*BTx1x+xp2'*BT2x'*BT1x];
P22=[0,0,0,vx1pr'*A1'*BTx2x+xp1'*BT1x'*BT2x];
end 
end
        
end




