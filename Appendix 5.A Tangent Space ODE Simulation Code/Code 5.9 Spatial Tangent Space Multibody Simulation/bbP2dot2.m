function [P21,P22]=bbP2dot2(i,j,a2pr,s1pr,s2pr,tn,q,x,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);
[r1,p1]=qPart(q,i);
BT1=BTran(p1,s1pr);
[xr1,xp1]=xPart(x,i);
BTx1=BTran(xp1,s1pr);

if j==0
P21=-[0,0,0,a2pr'*BTx1];
P22=zeros(1,7);    
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
A1=ATran(p1);
d12=r2+A2*s2pr-r1-A1*s1pr;
[xr2,xp2]=xPart(x,j);
BT2s2=BTran(p2,s2pr);
BT2a2=BTran(p2,a2pr);
BTx2s2=BTran(xp2,s2pr);
BTx2a2=BTran(xp2,a2pr);
c=xp2'*BT2a2'*BT2s2+d12'*BTx2a2+a2pr'*A2'*BTx2s2+(xr2'+xp2'*BT2s2')*BT2a2;
P21=-[xp2'*BT2a2',a2pr'*A2'*BT1+xp2'*BT2a2'*BT1];
P22=[xp2'*BT2a2',c-(xr1'+xp1'*BT1')*BT2a2];
end   
        
end






