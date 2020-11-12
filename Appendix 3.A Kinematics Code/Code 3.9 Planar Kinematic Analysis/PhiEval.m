 function Phi=PhiEval(q,PJDT,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

Phi=zeros(ngc,1);  %Set system constraint vector to zero

I2=eye(2);
P=[0,-1;1,0];
k=1;
m=0;
while k<=nh+nd  %Cycle through nh holonomic and nd driver constraints
    
if PJDT(1,k)==1     %Revolute constraint
[i,j,s1pr,s2pr]=RevPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1  %If body j is not ground
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
PhiRev=r2+A2*s2pr-r1-A1*s1pr;  %Evaluate revolute constraint
Phi=Add(Phi,PhiRev,m,0);  %Add revolute contribution to Phi
m=m+2;
end

if PJDT(1,k)==2     %Translational constrint
[i,j,s1pr,s2pr,v1pr,v2pr]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1  %If body j is not ground
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
PhiTran=-[v1pr'*P*A1'*d12;v1pr'*P*A1'*A2*v2pr];
Phi=Add(Phi,PhiTran,m,0);
m=m+2;
end

if PJDT(1,k)==3     %Distance constraint
[i,j,s1pr,s2pr,d]=DistPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1  %If body j is not ground
[r2,ph2]=qPart(q,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
PhiDist=(d12'*d12-d^2)/2;
Phi=Add(Phi,PhiDist,m,0);
m=m+1;
end

if PJDT(1,k)==4     %Relative rotation driver
[i,j]=RelRotPart(k,PJDT);
[r1,ph1]=qPart(q,i);
r2=[0;0];
ph2=0;
if j>=1
[r2,ph2]=qPart(q,j);
end
PhiRelRot=ph2-ph1;
Phi=Add(Phi,PhiRelRot,m,0);
m=m+1;
end

k=k+1;

end




