function [P21,P22]=bbP2dot2(i,j,a2pr,s1pr,s2pr,q,qd,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

I3=eye(3);
[r1,p1]=qPart(q,i);
A1=ATran(p1);
BT1s1=BTran(p1,s1pr);
[xr1,xp1]=qPart(qd,i);
BTx1s1=BTran(xp1,s1pr);

if j==0
P21=-[0,0,0,a2pr'*BTx1s1];
P22=zeros(1,7);    
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
d12=r2+A2*s2pr-r1-A1*s1pr;
[xr2,xp2]=qPart(qd,j);
BT2s2=BTran(p2,s2pr);
BT2a2=BTran(p2,a2pr);
BTx2s2=BTran(xp2,s2pr);
BTx2a2=BTran(xp2,a2pr);
c=xp2'*BT2a2'*BT2s2+d12'*BTx2a2+a2pr'*A2'*BTx2s2+(xr2'+xp2'*BT2s2')*BT2a2;
P21=-[xp2'*BT2a2',a2pr'*A2'*BTx1s1+xp2'*BT2a2'*BT1s1];
P22=[xp2'*BT2a2',c-(xr1'+xp1'*BT1s1')*BT2a2];
end   
        
end






