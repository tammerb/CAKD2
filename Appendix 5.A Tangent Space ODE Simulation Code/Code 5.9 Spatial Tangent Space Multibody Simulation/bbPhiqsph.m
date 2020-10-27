function [Phiq1,Phiq2]=bbPhiqsph(i,j,s1pr,s2pr,tn,q,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);

[r1,p1]=qPart(q,i);

if j==0
Phiq1=-[I3,BTran(p1,s1pr)];
Phiq2=zeros(3,7);
end

if j>=1
[r2,p2]=qPart(q,j);
Phiq1=-[I3,BTran(p1,s1pr)];
Phiq2=[I3,BTran(p2,s2pr)];
end

end





