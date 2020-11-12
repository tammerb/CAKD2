function [P31,P32]=bbP3dot1(i,j,a1pr,a2pr,tn,q,qd,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[r1,p1]=qPart(q,i);
[r1d,p1d]=qPart(qd,i);
BT1=BTran(p1,a1pr);
BT1d=BTran(p1d,a1pr);

if j==0
P31=zeros(1,7);
P32=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
BT2=BTran(p2,a2pr);
BT2d=BTran(p2d,a2pr);
P31=[0,0,0,2*p2d'*BT2'*BT1d+p2d'*BT2d'*BT1];
P32=[0,0,0,2*p1d'*BT1'*BT2d+p1d'*BT1d'*BT2];
end   
        
end




