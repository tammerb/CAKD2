function [Phiq1,Phiq2]=bbPhiqdot1(i,j,a1pr,a2pr,tn,q,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
A1=ATran(p1);
BT1=BTran(p1,a1pr);

if j==0
Phiq1=[0,0,0,a2pr'*BT1];
Phiq2=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
BT2=BTran(p2,a2pr);
Phiq1=[0,0,0,a2pr'*A2'*BT1];
Phiq2=[0,0,0,a1pr'*A1'*BT2];
end   
        
end
