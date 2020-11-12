function Phiq=PhiqEval(tn,q,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

Phiq=zeros(nc,ngc);

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
Phiksq1=-[I2,P*A1*s1pr];
Phiq=Add(Phiq,Phiksq1,m,3*(i-1));
if j>=1
Phiksq2=[I2,P*A2*s2pr];
Phiq=Add(Phiq,Phiksq2,m,3*(j-1));
end
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
Phiksq1=[v1pr'*P*A1',-v1pr'*s1pr-v1pr'*A1'*d12;0,0,-v1pr'*A1'*A2*v2pr];
Phiq=Add(Phiq,Phiksq1,m,3*(i-1));
if j>=1
Phiksq2=[-v1pr'*P*A1',v1pr'*A1'*A2*s2pr;0,0,v1pr'*A1'*A2*v2pr];
Phiq=Add(Phiq,Phiksq2,m,3*(j-1));
end
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
Phiksq1=-[d12',d12'*P*A1*s1pr];
Phiq=Add(Phiq,Phiksq1,m,3*(i-1));
if j>=1
Phiksq2=[d12',d12'*P*A2*s2pr];
Phiq=Add(Phiq,Phiksq2,m,3*(j-1));
end
m=m+1;
end

k=k+1;

end


