function [P31,P32]=bbP3dot2(i,j,a2pr,s1pr,s2pr,tn,q,qd,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);
[r1,p1]=qPart(q,i);
[r1d,p1d]=qPart(qd,i);
BT1=BTran(p1,s1pr);
BT1d=BTran(p1d,s1pr);

if j==0
P31=zeros(1,7);
P32=zeros(1,7);    
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
BT2s2=BTran(p2,s2pr);
BT2a2=BTran(p2,a2pr);
BT2s2d=BTran(p2d,s2pr);
BT2a2d=BTran(p2d,a2pr);
d=p2d'*BT2a2d'*BT2s2+p2d'*BT2s2d'*BT2a2+2*r2d'*BT2a2d+...
    2*p2d'*BT2s2'*BT2a2d+2*p2d'*BT2a2'*BT2s2d;
P31=-[p2d'*BT2a2',2*p2d'*BT2a2'*BT1d+p2d'*BT2a2d'*BT1];
P32=[p2d'*BT2a2',-p1d'*BT1d'*BT2a2-2*(r1d'+p1d'*BT1')*BT2a2d+d];
end   
        
end

