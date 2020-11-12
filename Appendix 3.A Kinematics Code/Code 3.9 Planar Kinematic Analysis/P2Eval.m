function P2=P2Eval(q,qd,PJDT,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

P2=zeros(ngc,ngc);

I2=eye(2);
Z2=zeros(2,2);
P=[0,-1;1,0];
k=1;        %Constraint Number
m=0;        %Row Index in P2 
while k<=nh+nd
    
if PJDT(1,k)==1     %Revolute
[i,j,s1pr,s2pr]=RevPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
if j==0
P21=[Z2,ph1d*A1*s1pr];   
P2=Add(P2,P21,m,3*(i-1));        
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
P21=[Z2,ph1d*A1*s1pr];
P2=Add(P2,P21,m,3*(i-1));
P22=[Z2,-ph2d*A2*s2pr];
P2=Add(P2,P22,m,3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==2     %Translational
[i,j,s1pr,s2pr,v1pr,v2pr]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
if j==0
r2=[0;0];
ph2=0;
r2d=[0;0];
ph2d=0;
A2=I2;
d12=r2+A2*s2pr-r1-A1*s1pr;
P21=[ph1d*v1pr'*A1',v1pr'*A1'*(r1d-r2d+ph1d*P*d12-ph2d*P*A2*s2pr+...
    ph1d*v1pr'*P*s1pr);0,0,(ph1d-ph2d)*v1pr'*A1'*P*A2*v2pr];
P2=Add(P2,P21,m,3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
P21=[ph1d*v1pr'*A1',v1pr'*A1'*(r1d-r2d+ph1d*P*d12-ph2d*P*A2*s2pr+...
    ph1d*v1pr'*P*s1pr);0,0,(ph1d-ph2d)*v1pr'*A1'*P*A2*v2pr];
P2=Add(P2,P21,m,3*(i-1));
P22=[-ph1d*v1pr'*A1',(ph2d-ph1d)*v1pr'*A1'*P*A2*s2pr;...
    0,0,(ph2d-ph1d)*v1pr'*A1'*P*A2*v2pr];
P2=Add(P2,P22,m,3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==3     %Distance
[i,j,s1pr,s2pr,d]=DistPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
if j==0
r2=[0;0];
ph2=0;
r2d=[0;0];
ph2d=0;
d12=s2pr-r1-A1*s1pr;
a1=r1d+ph1d*P*s1pr;
P21=[a1',a1'*P*A1*s1pr+ph1d*d12'*A1*s1pr];
P2=Add(P2,P21,m,3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
a1=r1d+ph1d*P*A1*s1pr;
a2=r2d+ph2d*P*A2*s2pr;
P21=[(a1-a2)',(a1-a2)'*P*A1*s1pr+ph1d*d12'*A1*s1pr];
P2=Add(P2,P21,m,3*(i-1));
P22=[(a2-a1)',(a2-a1)'*P*A2*s2pr+ph2d*d12'*A2*s2pr];
P2=Add(P2,P22,m,3*(j-1));
end
m=m+1;
end

k=k+1;
end
end




