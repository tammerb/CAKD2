function [P411,P412,P422]=bbP4dot1(i,j,a1pr,a2pr,tn,q,etak,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
BT1=BTran(p1,a1pr);
A1=ATran(p1);

if j==0
P411=etak*[zeros(3,7);zeros(4,3),KEval(a1pr,a2pr)];
P412=zeros(7,7);
P422=zeros(7,7);
end

if j>=1
[r2,p2]=qPart(q,j);
BT2=BTran(p2,a2pr);
A2=ATran(p2);

P411=etak*[zeros(3,7);zeros(4,3),KEval(a1pr,A2*a2pr)];
P412=etak*[zeros(3,7);zeros(4,3),BT1'*BT2];
P422=etak*[zeros(3,7);zeros(4,3),KEval(a2pr,A1*a1pr)];
end
         
end


