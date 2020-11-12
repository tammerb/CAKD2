function [P31,P32]=bbP3dist(i,j,s1pr,s2pr,d,tn,q,qd,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);

[r1,p1]=qPart(q,i);
[r1d,p1d]=qPart(qd,i);
BT1=BTran(p1,s1pr);
BT1d=BTran(p1d,s1pr);
a1bar=r1d+BT1*p1d;

if j==0
P31=[p1d'*BT1d',2*(a1bar)'*BT1d+p1d'*BT1d'*BT1];
P32=zeros(1,7);
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
BT2=BTran(p2,s2pr);    
BT2d=BTran(p2d,s2pr);
a2bar=r2d+BT2*p2d;
P31=[p1d'*BT1d'-p2d'*BT2d',2*(a1bar-a2bar)'*BT1d+...
    (p1d'*BT1d'-p2d'*BT2d')*BT1];
P32=[p2d'*BT2d'-p1d'*BT1d',2*(a2bar-a1bar)'*BT2d+...
    (p2d'*BT2d'-p1d'*BT1d')*BT2];
end
end

