function P2=P2Eval(tn,q,x,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

P2=zeros(nc,ngc);

I2=eye(2);
Z2=zeros(2,2);
P=[0,-1;1,0];
k=1;        %Constraint Number
m=0;        %Row Index in P2 
while k<=nh
    
if PJDT(1,k)==1     %Revolute
[i,j,s1pr,s2pr]=RevPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[xr1,xph1]=xPart(x,i);
A1=ATran(ph1);
if j==0
P21=[Z2,xph1*A1*s1pr];   
P2=Add(P2,P21,m,3*(i-1));        
end
if j>=1
[r2,ph2]=qPart(q,j);
[xr2,xph2]=xPart(x,j);
A2=ATran(ph2);
P21=[Z2,xph1*A1*s1pr];
P2=Add(P2,P21,m,3*(i-1));
P22=[Z2,-xph2*A2*s2pr];
P2=Add(P2,P22,m,3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==2     %Translational
[i,j,s1pr,s2pr,v1pr,v2pr]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[xr1,xph1]=xPart(x,i);
A1=ATran(ph1);
if j==0
r2=[0;0];
ph2=0;
xr2=[0;0];
xph2=0;
A2=I2;
d12=r2+A2*s2pr-r1-A1*s1pr;
P21=[xph1*v1pr'*A1',v1pr'*A1'*(xr1-xr2+xph1*P*d12-xph2*P*A2*s2pr+...
    xph1*v1pr'*P*s1pr);0,0,(xph1-xph2)*v1pr'*A1'*P*A2*v2pr];
P2=Add(P2,P21,m,3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[xr2,xph2]=xPart(x,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
P21=[xph1*v1pr'*A1',v1pr'*A1'*(xr1-xr2+xph1*P*d12-xph2*P*A2*s2pr+...
    xph1*v1pr'*P*s1pr);0,0,(xph1-xph2)*v1pr'*A1'*P*A2*v2pr];
P2=Add(P2,P21,m,3*(i-1));
P22=[-xph1*v1pr'*A1',(xph2-xph1)*v1pr'*A1'*P*A2*s2pr;...
    0,0,(xph2-xph1)*v1pr'*A1'*P*A2*v2pr];
P2=Add(P2,P22,m,3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==3     %Distance
[i,j,s1pr,s2pr,d]=DistPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[xr1,xph1]=xPart(x,i);
A1=ATran(ph1);
if j==0
r2=[0;0];
ph2=0;
xr2=[0;0];
xph2=0;
d12=s2pr-r1-A1*s1pr;
a1=xr1+xph1*P*s1pr;
P21=[a1',a1'*P*A1*s1pr+xph1*d12'*A1*s1pr];
P2=Add(P2,P21,m,3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[xr2,xph2]=xPart(x,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
a1=xr1+xph1*P*A1*s1pr;
a2=xr2+xph2*P*A2*s2pr;
P21=[(a1-a2)',(a1-a2)'*P*A1*s1pr+xph1*d12'*A1*s1pr];
P2=Add(P2,P21,m,3*(i-1));
P22=[(a2-a1)',(a2-a1)'*P*A2*s2pr+xph2*d12'*A2*s2pr];
P2=Add(P2,P22,m,3*(j-1));
end
m=m+1;
end

k=k+1;
end
end




