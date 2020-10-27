function P3=P3Eval(tn,q,qd,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

P3=zeros(nc,ngc);

I2=eye(2);
Z2=zeros(2,2);
P=[0,-1;1,0];
k=1;        %Constraint Number
m=0;        %Row Index in P3 
while k<=nh
    
if PJDT(1,k)==1     %Revolute
[i,j,s1pr,s2pr]=RevPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
if j==0
P31=[Z2,(ph1d^2)*P*A1*s1pr];
P3=Add(P3,P31,m,3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
P31=[Z2,(ph1d^2)*P*A1*s1pr];
P3=Add(P3,P31,m,3*(i-1));
P32=[Z2,-(ph2d^2)*A2*s2pr];
P3=Add(P3,P32,m,3*(j-1));
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
P31=[-(ph1d^2)*v1pr'*A1'*P,(ph1d^2)*v1pr'*s1pr+...
    2*ph1d*v1pr'*A1'*P*(r2d-r1d)+(ph1d^2)*v1pr'*A1'*d12+...
    ((ph2d^2)-2*ph1d*ph2d)*v1pr'*A1'*A2*s2pr;...
    0,0,((ph1d^2)-2*ph1d*ph2d+(ph2d^2))*v1pr'*A1'*A2*v2pr];
P3=Add(P3,P31,m,3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
P31=[-(ph1d^2)*v1pr'*A1'*P,(ph1d^2)*v1pr'*s1pr+...
    2*ph1d*v1pr'*A1'*P*(r2d-r1d)+(ph1d^2)*v1pr'*A1'*d12+...
    ((ph2d^2)-2*ph1d*ph2d)*v1pr'*A1'*A2*s2pr;...
    0,0,((ph1d^2)-2*ph1d*ph2d+(ph2d^2))*v1pr'*A1'*A2*v2pr];
P3=Add(P3,P31,m,3*(i-1));
P32=[(ph1d^2)*v1pr'*A1'*P,-...
    ((ph1d^2)-2*ph1d*ph2d+(ph2d^2))*v1pr'*A1'*A2*s2pr;...
    0,0,-((ph1d^2)-2*ph1d*ph2d+(ph2d^2))*v1pr'*A1'*A2*v2pr];
P3=Add(P3,P32,m,3*(j-1));
end
m=m+2;
end

if PJDT(1,k)==3     %Distance
[i,j,s1pr,s2pr,d]=DistPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
ab1=r1d+ph1d*P*A1*s1pr;
if j==0
r2=[0;0];
ph2=0;
r2d=[0;0];
ph2d=0;    
A2=I2;
d12=r2+A2*s2pr-r1-A1*s1pr;    
ab2=r2d+ph2d*P*A2*s2pr; 
P31=-s1pr'*A1'*[(ph1d^2)*I2,2*ph1d*ab1+(ph1d^2)*P*d12]+...
    2*[0,0,ph1d*ab2'*A1*s1pr]+(ph2d^2)*s2pr'*A2'*[I2,P*A1*s1pr];
P3=Add(P3,P31,m,3*(i-1));   
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;    
ab2=r2d+ph2d*P*A2*s2pr; 
P31=-s1pr'*A1'*[(ph1d^2)*I2,2*ph1d*ab1+(ph1d^2)*P*d12]+...
    2*[0,0,ph1d*ab2'*A1*s1pr]+(ph2d^2)*s2pr'*A2'*[I2,P*A1*s1pr];
P3=Add(P3,P31,m,3*(i-1));
P32=-s2pr'*A2'*[(ph2d^2)*I2,2*ph2d*ab2-(ph2d^2)*P*d12]+...
    2*[0,0,ph2d*ab1'*A2*s2pr]+(ph1d^2)*s1pr'*A1'*[I2,P*A2*s2pr];
P3=Add(P3,P32,m,3*(j-1));
end
m=m+1;
end

k=k+1;
end

end





