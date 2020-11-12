function [P411,P412,P422]=bbP4dist(i,j,s1pr,s2pr,d,tn,q,etak,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);

[r1,p1]=qPart(q,i);
BT1=BTran(p1,s1pr);
A1=ATran(p1);

if j==0
d12=s2pr-r1-A1*s1pr;
P411=etak*[I3,BT1;BT1',BT1'*BT1-KEval(s1pr,d12)];
P412=zeros(7,7);
P422=zeros(7,7);
end

if j>=1
[r2,p2]=qPart(q,j);
BT2=BTran(p2,s2pr);
A2=ATran(p2);
d12=r2+A2*s2pr-r1-A1*s1pr;
P411=etak*[I3,BT1;BT1',BT1'*BT1-KEval(s1pr,d12)];
P412=-etak*[I3,BT2;BT1',BT1'*BT2];
P422=etak*[I3,BT2;BT2',BT2'*BT2+KEval(s2pr,d12)];
end
         
end

