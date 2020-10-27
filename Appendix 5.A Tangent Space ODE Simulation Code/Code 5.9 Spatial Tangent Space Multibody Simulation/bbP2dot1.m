function [P21,P22]=bbP2dot1(i,j,a1pr,a2pr,tn,q,x,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
[xr1,xp1]=xPart(x,i);
A1=ATran(p1);
BT1=BTran(p1,a1pr);
BTx1=BTran(xp1,a1pr);

if j==0
P21=[0,0,0,a2pr'*BTx1];
P22=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
[xr2,xp2]=xPart(x,j);
A2=ATran(p2);
BT2=BTran(p2,a2pr);
BTx2=BTran(xp2,a2pr);
P21=[0,0,0,a2pr'*A2'*BTx1+xp2'*BT2'*BT1];
P22=[0,0,0,a1pr'*A1'*BTx2+xp1'*BT1'*BT2];
end   
        
end


