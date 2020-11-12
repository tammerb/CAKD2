function P4=P4Eval(tn,q,eta,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

P4=zeros(ngc,ngc);

I2=eye(2);
Z1=zeros(2,1);
Z2=zeros(2,2);
P=[0,-1;1,0];
k=1;    %Joint number
m=0;    %Address in vector eta

while k<=nh
    
if PJDT(1,k)==1     %Revolute
[i,j,s1pr,s2pr]=RevPart(k,PJDT);
[r1,ph1]=qPart(q,i);
A1=ATran(ph1);
etak=[eta(m+1);eta(m+1)];
if j==0
P411=[Z2,Z1;0,0,-s1pr'*A1'*etak];   
P4=Add(P4,P411,3*(i-1),3*(i-1));       
end
if j>=1
[r2,ph2]=qPart(q,j);
A2=ATran(ph2);
P411=[Z2,Z1;0,0,-s1pr'*A1'*etak];
P4=Add(P4,P411,3*(i-1),3*(i-1));
P422=[Z2,Z1;0,0,-s2pr'*A2'*etak];
P4=Add(P4,P422,3*(j-1),3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==2     %Translational
[i,j,s1pr,s2pr,v1pr,v2pr]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
A1=ATran(ph1);
eta1=eta(m+1);
eta2=eta(m+2);
if j==0
r2=[0;0];
A2=I2;
d12=r2+A2*s2pr-r1-A1*s1pr;
P411=[Z2,eta1*A1*v1pr;eta1*v1pr'*A1',eta1*v1pr'*P*s1pr-...
    eta1*d12'*P*A1*v1pr+eta2*v2pr'*A2'*P*A1*v1pr];
P4=Add(P4,P411,3*(i-1),3*(i-1));
end    
if j>=1
[r2,ph2]=qPart(q,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
P411=[Z2,eta1*A1*v1pr;eta1*v1pr'*A1',eta1*v1pr'*P*s1pr-...
    eta1*d12'*P*A1*v1pr+eta2*v2pr'*A2'*P*A1*v1pr];
P4=Add(P4,P411,3*(i-1),3*(i-1));
P412=[Z2,Z1;-eta1*v1pr'*A1',-v1pr'*A1'*P*A2*(eta1*s1pr-eta2*v2pr)];
P421=P412';
P422=[Z2,Z1;0,0,v1pr'*A1'*P*A2*(eta1*s2pr+eta2*v2pr)];
P4=Add(P4,P412,3*(i-1),3*(j-1));
P4=Add(P4,P421,3*(j-1),3*(i-1));
P4=Add(P4,P422,3*(j-1),3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==3     %Distance
[i,j,s1pr,s2pr,d]=DistPart(k,PJDT);
[r1,ph1]=qPart(q,i);
A1=ATran(ph1);
etak=eta(m+1);
if j==0
d12=-r1-A1*s1pr;    
P411=etak*[I2,P*A1*s1pr;-s1pr'*A1'*P,s1pr'*s1pr+s1pr'*A1'*d12];
P4=Add(P4,P411,3*(i-1),3*(i-1));   
end
if j>=1
[r2,ph2]=qPart(q,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
P411=etak*[I2,P*A1*s1pr;-s1pr'*A1'*P,s1pr'*s1pr+s1pr'*A1'*d12];
P4=Add(P4,P411,3*(i-1),3*(i-1));
P412=etak*[-I2,-P*A2*s2pr;s1pr'*A1'*P,-s1pr'*A1'*A2*s2pr];
P421=P412;
P422=etak*[I2,P*A2*s2pr;-s2pr'*A2'*P,s2pr'*s2pr-s2pr'*A2'*d12];   
P4=Add(P4,P412,3*(i-1),3*(j-1));
P4=Add(P4,P421,3*(j-1),3*(i-1));
P4=Add(P4,P422,3*(j-1),3*(j-1));
end
m=m+1;
end

k=k+1;
end
end




