function [P21,P22]=bbP2sph(i,j,s1pr,s2pr,tn,q,x,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

[xr1,xp1]=xPart(x,i);
BT1x=BTran(xp1,s1pr);

if j==0
P21=[zeros(3,3),-BT1x];
P22=zeros(3,7);        
end

if j>=1
[xr2,xp2]=xPart(x,j);
BT2x=BTran(xp2,s2pr);
P21=[zeros(3,3),-BT1x];   
P22=[zeros(3,3),BT2x];
end

end







