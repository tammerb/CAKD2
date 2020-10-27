function Phi=PhiEval(tn,q,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

Phi=zeros(nc,1);

I2=eye(2);
P=[0,-1;1,0];
k=1;
m=0;
while k<=nh
    
if PJDT(1,k)==1     %Revolute
[i,j,s1pr,s2pr]=RevPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
PhiRev=r2+A2*s2pr-r1-A1*s1pr;
Phi=Add(Phi,PhiRev,m,0);
m=m+2;
end

if PJDT(1,k)==2     %Translational
[i,j,s1pr,s2pr,v1pr,v2pr]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
PhiTran=-[v1pr'*P*A1'*d12;v1pr'*P*A1'*A2*v2pr];
Phi=Add(Phi,PhiTran,m,0);
m=m+2;
end

if PJDT(1,k)==3     %Distance
[i,j,s1pr,s2pr,d]=DistPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
PhiDist=(d12'*d12-d^2)/2;
Phi=Add(Phi,PhiDist,m,0);
m=m+1;
end

k=k+1;

end




