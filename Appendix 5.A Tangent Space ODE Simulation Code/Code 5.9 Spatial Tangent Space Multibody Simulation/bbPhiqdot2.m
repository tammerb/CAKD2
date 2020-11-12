function [Phiq1,Phiq2]=bbPhiqdot2(i,j,a2pr,s1pr,s2pr,tn,q,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);
[r1,p1]=qPart(q,i);
A1=ATran(p1);
BT1=BTran(p1,s1pr);

if j==0
d12=s2pr-r1-A1*s1pr;
Phiq1=-a2pr'*[I3,BT1];
Phiq2=zeros(1,7);    
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
BT2s2=BTran(p2,s2pr);
BT2a2=BTran(p2,a2pr);
d12=r2+A2*s2pr-r1-A1*s1pr;
Phiq1=-a2pr'*A2'*[I3,BT1];
Phiq2=a2pr'*A2'*[I3,BT2s2]+d12'*[zeros(3,3),BT2a2];
end   
        
end




