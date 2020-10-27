function [P21,P22]=bbP2dist(i,j,s1pr,s2pr,d,tn,q,x,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
[xr1,xp1]=xPart(x,i);
A1=ATran(p1);
BT1=BTran(p1,s1pr);
BT1x=BTran(xp1,s1pr);
a1=(xr1'+xp1'*BT1');

if j==0
d12=s2pr-r1-A1*s1pr;
P21=[a1,a1*BT1-d12'*BT1x];
P22=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
[xr2,xp2]=xPart(x,j);
A2=ATran(p2);
BT2=BTran(p2,s2pr);
BT2x=BTran(xp2,s2pr);
a2=(xr2'+xp2'*BT2');
d12=r2+A2*s2pr-r1-A1*s1pr;
P21=[(a1-a2),(a1-a2)*BT1-d12'*BT1x];
P22=[(a2-a1),(a2-a1)*BT2+d12'*BT2x];
end
         
end




