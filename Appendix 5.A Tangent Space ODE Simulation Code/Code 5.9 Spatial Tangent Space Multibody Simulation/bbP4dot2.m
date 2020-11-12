function [P411,P412,P422]=bbP4dot2(i,j,a2pr,s1pr,s2pr,tn,q,etak,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
BT1=BTran(p1,s1pr);
A1=ATran(p1);

if j==0
P411=-etak*[zeros(3,7);zeros(4,3),KEval(s1pr,a2pr)];
P412=zeros(7,7);
P422=zeros(7,7);
end

if j>=1
[r2,p2]=qPart(q,j);
A2=ATran(p2);
BT2a2=BTran(p2,a2pr);
BT2s2=BTran(p2,s2pr);
d12=r2+A2*s2pr-r1-A1*s1pr;
e=KEval(s2pr,A2*a2pr)+KEval(a2pr,d12)+BT2s2'*BT2a2+BT2a2'*BT2s2;
P411=-etak*[zeros(3,7);zeros(4,3),KEval(s1pr,A2*a2pr)];
P412=-etak*[zeros(3,3),BT2a2;zeros(4,3),BT1'*BT2a2];
P422=etak*[zeros(3,3),BT2a2;BT2a2',e];
end
         
end




