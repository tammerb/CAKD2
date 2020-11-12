function [QAsq,QAsqd]=QAsqqdEval(tn,q,qd,Lam,...
    PMDT,PJDT,PTSDAT,PRSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA,NRSDA]=...
    parPart(par);

uy=[0;1];

QAsq=zeros(ngc,ngc);
QAsqd=zeros(ngc,ngc);

P=[0,-1;1,0];
I2=eye(2);
Z2=zeros(2);

%Account for TSDA forces
T=1;
while T<=NTSDA
    
%Evaluate QAsq
[i,j,s1pr,s2pr,K,C,el0,F]=PTSDATPart(PTSDAT,T);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
r2=[0;0];
ph2=0;
r2d=[0;0];
ph2d=0;
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
end
A1=ATran(ph1);
A2=ATran(ph2);
d12=r2+A2*s2pr-r1-A1*s1pr;
el=sqrt(d12'*d12);
a=r2d+ph2d*P*A2*s2pr-r1d-ph1d*P*A1*s1pr;
eld=(1/el)*d12'*a;

%K=K*(1-sign(q(1)))/2;     %Unilateral spring constant
f=K*(el-el0)+C*eld+F;   %User insert F(el,eld) if needed,
elsq1=-(1/el)*d12'*[I2,P*A1*s1pr];
eldsq1=-(1/el^2)*d12'*a*elsq1-(1/el)*a'*[I2,P*A1*s1pr]+...
    (1/el)*d12'*[Z2,ph1d*A1*s1pr];

%K=K*(1-sign(q(1)))/2;     %Unilateral spring constant
b1=(1/el)*(K*elsq1+C*eldsq1)-(f/el^2)*elsq1;
Q1sq1=[d12;-s1pr'*A1'*P*d12]*b1+(f/el)*[[-I2,-P*A1*s1pr];...
    [-s1pr'*A1'*P,-s1pr'*s1pr-s1pr'*A1'*d12]];
QAsq=Add(QAsq,Q1sq1,3*(i-1),3*(i-1));

if j>=1
elsq2=(1/el)*d12'*[I2,P*A2*s2pr];
eldsq2=-(1/el^2)*d12'*a*elsq2+(1/el)*a'*[I2,P*A2*s2pr]-...
    (1/el)*d12'*[Z2,ph2d*A2*s2pr];

%K=K*(1-sign(q(1)))/2;     %Unilateral spring constant
b2=(1/el)*(K*elsq2+C*eldsq2)-(f/el^2)*elsq2;    
Q1sq2=[d12;-s1pr'*A1'*P*d12]*b2+(f/el)*[[I2,P*A2*s2pr];...
    [-s1pr'*A1'*P*[I2,P*A2*s2pr]]];
QAsq=Add(QAsq,Q1sq2,3*(i-1),3*(j-1));
Q2sq1=-[d12;-s2pr'*A2'*P*d12]*b1-(f/el)*[[-I2,-P*A1*s1pr];...
    [s2pr'*A2'*P*[I2,P*A1*s1pr]]];
QAsq=Add(QAsq,Q2sq1,3*(j-1),3*(i-1));
Q2sq2=-[d12;-s2pr'*A2'*P*d12]*b2-(f/el)*[[I2,P*A2*s2pr];...
    [-s2pr'*A2'*P,-s2pr'*s2pr+s2pr'*A2'*d12]];
QAsq=Add(QAsq,Q2sq2,3*(j-1),3*(j-1));
end

%Evaluate QAsqd
eldsq1d=(1/el)*d12'*[-I2,-P*A1*s1pr];
Q1sq1d=(C/el)*[d12;-s1pr'*A1'*P*d12]*eldsq1d;
QAsqd=Add(QAsqd,Q1sq1d,3*(i-1),3*(i-1));

if j>=1
eldsq2d=(1/el)*d12'*[I2,-P*A2*s2pr];    
Q1sq2d=(C/el)*[d12;s1pr'*A1'*P*d12]*eldsq2d;
QAsqd=Add(QAsqd,Q1sq2d,3*(i-1),3*(j-1));
Q2sq1d=-(C/el)*[d12;-s2pr'*A2'*P*d12]*eldsq1d;
QAsqd=Add(QAsqd,Q2sq1d,3*(j-1),3*(i-1));
Q2sq2d=-(C/el)*[d12;-s2pr'*A2'*P*d12]*eldsq2d;
QAsqd=Add(QAsqd,Q2sq2d,3*(j-1),3*(j-1));
end   

T=T+1;
end

%Account for RSDA forces

m=1;
while m<=NRSDA
[i,j,K,C,phi0,T]=PRSDATPart(PRSDAT,m);

%Evaluate QAsq
QA1sq1=[zeros(2,3);0,0,-K];
QAsq=Add(QAsq,QA1sq1,3*(i-1),3*(i-1));

if j>=1
QA1sq2=[zeros(2,3);0,0,K];
QA2sq1=[zeros(2,3);0,0,K];
QA2sq2=[zeros(2,3);0,0,-K];
QAsq=Add(QAsq,QA1sq2,3*(i-1),3*(j-1));
QAsq=Add(QAsq,QA2sq1,3*(j-1),3*(i-1));
QAsq=Add(QAsq,QA2sq2,3*(j-1),3*(j-1));
end

%Evaluate QAsqd
QA1sq1d=[zeros(2,3);0,0,-C];
QAsqd=Add(QAsqd,QA1sq1d,3*(i-1),3*(i-1));

if j>=1
QA1sq2d=[zeros(2,3);0,0,C];
QA2sq1d=[zeros(2,3);0,0,C];
QA2sq2d=[zeros(2,3);0,0,-C];
QAsqd=Add(QAsqd,QA1sq2d,3*(i-1),3*(j-1));
QAsqd=Add(QAsqd,QA2sq1d,3*(j-1),3*(i-1));
QAsqd=Add(QAsqd,QA2sq2d,3*(j-1),3*(j-1));
end
m=m+1;
end

%Account for terms due to friction in Revolute and Translational joints
k=1;
while k<=nh

if PJDT(1,k)==1     %Revolute joint
[i,j,s1pr,s2pr,R,mus,mud,ms,nm]=RevPart(k,PJDT);
[r1d,ph1d]=qPart(qd,i);
Lamk=[Lam(ms);Lam(ms+1)];
if j==0
[Sfr,Sfrpr]=SfrSfrpr(R*(-ph1d),mus,mud,par);
QAfk1sq1d=(R^2)*norm(Lamk)*[zeros(2,3);0,0,-Sfrpr];
QAsqd=Add(QAsqd,QAfk1sq1d,3*(i-1),3*(i-1));    
end
if j>=1
[r2d,ph2d]=qPart(qd,j);
[Sfr,Sfrpr]=SfrSfrpr(R*(ph2d-ph1d),mus,mud,par);
QAfk1sq1d=(R^2)*norm(Lamk)*[zeros(2,3);0,0,-Sfrpr];
QAsqd=Add(QAsqd,QAfk1sq1d,3*(i-1),3*(i-1));
QAfk1sq2d=-QAfk1sq1d;
QAsqd=Add(QAsqd,QAfk1sq2d,3*(i-1),3*(j-1));
QAfk2sq1d=-QAfk1sq1d;
QAsqd=Add(QAsqd,QAfk2sq1d,3*(j-1),3*(i-1));
QAfk2sq2d=-QAfk1sq2d;
QAsqd=Add(QAsqd,QAfk2sq2d,3*(j-1),3*(j-1));
end

end

if PJDT(1,k)==2     %Translational joint
[i,j,s1pr,s2pr,v1pr,v2pr,mus,mud,d,ms,nm]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
Lamk=[Lam(ms);Lam(ms+1)];
Lamk1=Lamk(1);
Lamk2=Lamk(2);
P=[0,-1;1,0];
if j==0
T1k=-s1pr'*v1pr*Lamk1+...
    [(s2pr'-r1')*A1*v1pr,v1pr'*A1'*v2pr]*Lamk;
f1k=(T1k/d)-d*Lamk1;
f2k=-T1k/d;
v12=v1pr'*A1'*(-r1d)/d-ph1d*v1pr'*P*s1pr/d;
[Sfr,Sfrpr]=SfrSfrpr(v12,mus,mud,par);
a1=[(s2pr'-r1')*P*A1*v1pr,-v1pr'*A1'*P*v2pr]*Lamk;
T1ksq1=[-v1pr'*A1',a1];
f1ksq1=(1/d)*T1ksq1;
f2ksq1=-f1ksq1;
[csf1k,dcsf1k]=csign(f1k,par);
[csf2k,dcsf2k]=csign(f2k,par);
absf1k=csf1k*f1k;
absf2k=csf2k*f2k;
absf1ksq1=(csf1k+f1k*dcsf1k)*f1ksq1;
absf2ksq1=(csf2k+f2k*dcsf2k)*f2ksq1;
v12sq1=[zeros(1,2),-(1/d)*v1pr'*A1'*P*(-r1d)];
v12sq1d=(1/d)*[-v1pr'*A1',-v1pr'*P*s1pr];
QAfk1sq1=(1/d)*Sfr*[A1*v1pr;v1pr'*P*s1pr]*(absf1ksq1+absf2ksq1)+...
    (1/d)*(absf1k+absf2k)*[A1*v1pr;v1pr'*P*s1pr]*Sfrpr*v12sq1+...
    (1/d)*(absf1k+absf2k)*Sfr*[zeros(2),P*A1*v1pr;zeros(1,3)];
QAsq=Add(QAsq,QAfk1sq1,3*(i-1),3*(i-1));
QAfk1sq1d=(1/d)*(absf1k+absf2k)*[A1*v1pr;v1pr'*P*s1pr]*Sfrpr*v12sq1d;
QAsqd=Add(QAsqd,QAfk1sq1d,3*(i-1),3*(i-1));
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
T1k=-s1pr'*v1pr*Lamk1+...
    [(r2'+s2pr'*A2'-r1')*A1*v1pr,v1pr'*A1'*A2*v2pr]*Lamk;
f1k=(T1k/d)-d*Lamk1;
f2k=-T1k/d;
v12=v1pr'*A1'*(r2d+ph2d*P*A2*s2pr-r1d)/d-ph1d*v1pr'*P*s1pr/d;
[Sfr,Sfrpr]=SfrSfrpr(v12,mus,mud,par);
a1=[(r2'+s2pr'*A2'-r1')*P*A1*v1pr,-v1pr'*A1'*P*A2*v2pr]*Lamk;
T1ksq1=[-v1pr'*A1',a1];
f1ksq1=(1/d)*T1ksq1;
f2ksq1=-f1ksq1;
[csf1k,dcsf1k]=csign(f1k,par);
[csf2k,dcsf2k]=csign(f2k,par);
absf1k=csf1k*f1k;
absf2k=csf2k*f2k;
absf1ksq1=(csf1k+f1k*dcsf1k)*f1ksq1;
absf2ksq1=(csf2k+f2k*dcsf2k)*f2ksq1;
v12sq1=[zeros(1,2),-(1/d)*v1pr'*A1'*P*(r2d+ph2d*P*A2*s2pr-r1d)];
v12sq1d=(1/d)*[-v1pr'*A1',-v1pr'*P*s1pr];
QAfk1sq1=(1/d)*Sfr*[A1*v1pr;v1pr'*P*s1pr]*(absf1ksq1+absf2ksq1)+...
    (1/d)*(absf1k+absf2k)*[A1*v1pr;v1pr'*P*s1pr]*Sfrpr*v12sq1+...
    (1/d)*(absf1k+absf2k)*Sfr*[zeros(2),P*A1*v1pr;zeros(1,3)];
QAsq=Add(QAsq,QAfk1sq1,3*(i-1),3*(i-1));
QAfk1sq1d=(1/d)*(absf1k+absf2k)*[A1*v1pr;v1pr'*P*s1pr]*Sfrpr*v12sq1d;
QAsqd=Add(QAsqd,QAfk1sq1d,3*(i-1),3*(i-1));
a2=[-s2pr'*A2'*P*A1*v1pr,v1pr'*A1'*P*A2*v2pr]*Lamk;
T1ksq2=[v1pr'*A1',a2];
f1ksq2=(1/d)*T1ksq2;
f2ksq2=-f1ksq2;
absf1ksq2=(csf1k+f1k*dcsf1k)*f1ksq2;
absf2ksq2=(csf2k+f2k*dcsf2k)*f2ksq2;
v12sq2=[zeros(1,2),-(1/d)*ph2d*v1pr'*A1'*A2*s2pr];
v12sq2d=(1/d)*[v1pr'*A1',v1pr'*A1'*P*A2*s2pr];
QAfk1sq2=(1/d)*Sfr*[A1*v1pr;v1pr'*P*s1pr]*(absf1ksq2+absf2ksq2)+...
    (1/d)*(absf1k+absf2k)*[A1*v1pr;v1pr'*P*s1pr]*Sfrpr*v12sq2;
QAsq=Add(QAsq,QAfk1sq2,3*(i-1),3*(j-1));
QAfk2sq1=(1/d)*Sfr*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*(absf1ksq1+absf2ksq1)+...
   (1/d)*(absf1k+absf2k)*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*Sfrpr*v12sq1+...
   (1/d)*(absf1k+absf2k)*Sfr*[zeros(2),-P*A1*v1pr;0,0,-v1pr'*A1'*A2*s2pr];
QAsq=Add(QAsq,QAfk2sq1,3*(j-1),3*(i-1));
QAfk2sq2=(1/d)*Sfr*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*(absf1ksq2+absf2ksq2)+...
   (1/d)*(absf1k+absf2k)*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*Sfrpr*v12sq2+...
   (1/d)*(absf1k+absf2k)*Sfr*[zeros(2,3);0,0,v1pr'*A1'*A2*s2pr];
QAsq=Add(QAsq,QAfk2sq2,3*(j-1),3*(j-1));
QAfk1sq2d=(1/d)*(absf1k+absf2k)*[A1*v1pr;v1pr'*P*s1pr]*Sfrpr*v12sq2d;
QAsqd=Add(QAsqd,QAfk1sq2d,3*(i-1),3*(j-1));
QAfk2sq1d=(1/d)*(absf1k+absf2k)*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*...
    Sfrpr*v12sq1d;
QAsqd=Add(QAsqd,QAfk2sq1d,3*(j-1),3*(i-1));
QAfk2sq2d=(1/d)*(absf1k+absf2k)*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*...
    Sfrpr*v12sq2d;
QAsqd=Add(QAsqd,QAfk2sq2d,3*(j-1),3*(j-1));   
end

end

k=k+1;
end

end



