function QAw=QAwEval(q,qd,Lam,PMDT,PJDT,PTSDAT,PRSDAT,par,w,N)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

uy=[0;1];

QA=zeros(ngc,1);

%Account for gravitational force in negative y direction
i=1;
while i<=nb
mi=PMDT(1,i);
QAgi=[-mi*g*uy;0];
QA=Add(QA,QAgi,3*(i-1),0);
i=i+1;
end

%Account for TSDA forces
P=[0,-1;1,0];
m=1;
while m<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=PTSDATPart(PTSDAT,m);
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
el=sqrt(d12'*d12+10^-4);
eld=(1/el)*d12'*(r2d+ph2d*P*A2*s2pr-r1d-ph1d*P*A1*s1pr);
%K=K*(1-sign(q(1)))/2;     %Unilateral spring constant
f=K*(el-el0)+C*eld+F;   %User insert F(el,eld) if needed, 
QA1=(f/el)*[d12;d12'*P*A1*s1pr];
QA=Add(QA,QA1,3*(i-1),0);
if j>=1
QA2=-(f/el)*[d12;d12'*P*A2*s2pr];
QA=Add(QA,QA2,3*(j-1),0);
end

m=m+1;
end

%Account for RSDA forces
m=1;
while m<=NRSDA
[i,j,K,C,phi0,T]=PRSDATPart(PRSDAT,m);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
ph2=0;
ph2d=0;
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
end
tau=K*(ph2-ph1-phi0)+C*(ph2d-ph1d)+T;   %User insert T(qd) if needed
QA1=[0;0;tau];
QA=Add(QA,QA1,3*(i-1),0);
if j>=1
QA2=-[0;0;tau];
QA=Add(QA,QA2,3*(j-1),0);
end

m=m+1;
end

%Account for friction in revolute and translational joints 
%and constraint force in distance constraints
k=1;
while k<=nh

if PJDT(1,k)==1     %Revolute joint;
[i,j,s1pr,s2pr,R,mus,mud,ms,nm]=RevPart(k,PJDT);
[r1d,ph1d]=qPart(qd,i);
Lamk=[Lam(ms);Lam(ms+1)];
if j==0
[Sfr,Sfrpr]=SfrSfrpr(R*(-ph1d),mus,mud,par);
Sfr=(w/N)*Sfr;
Sfrpr=(w/N)*Sfrpr;
QAfk1=[0;0;R*norm(Lamk)*Sfr];
QA=Add(QA,QAfk1,3*(i-1),0);
end
if j>=1
[r2d,ph2d]=qPart(qd,j);
[Sfr,Sfrpr]=SfrSfrpr(R*(ph2d-ph1d),mus,mud,par);
Sfr=(w/N)*Sfr;
Sfrpr=(w/N)*Sfrpr;
QAfk1=[0;0;R*norm(Lamk)*Sfr];
QA=Add(QA,QAfk1,3*(i-1),0);
QAfk2=-QAfk1;
QA=Add(QA,QAfk2,3*(j-1),0);
end
F(k)=norm(Lamk);
end

if PJDT(1,k)==2     %Translational joint
[i,j,s1pr,s2pr,v1pr,v2pr,mus,mud,d,ms,nm]=TranPart(k,PJDT);
[r1,ph1]=qPart(q,i);
[r1d,ph1d]=qPart(qd,i);
A1=ATran(ph1);
Lamk=[Lam(ms);Lam(ms+1)];
P=[0,-1;1,0];
if j==0
T1k=-s1pr'*v1pr*Lamk(1)+...
    [(s2pr'-r1')*A1*v1pr,v1pr'*A1'*v2pr]*Lamk;
f1k=(T1k/d)-d*Lamk(1);
f2k=-T1k/d;
v12=v1pr'*A1'*(-r1d)/d-ph1d*v1pr'*P*s1pr/d;
[Sfr,Sfrpr]=SfrSfrpr(v12,mus,mud,par);
Sfr=(w/N)*Sfr;
Sfrpr=(w/N)*Sfrpr;
[csf1k,dcsf1k]=csign(f1k,par);
[csf2k,dcsf2k]=csign(f2k,par);
QAfk1=(1/d)*(f1k*csf1k+f2k*csf2k)*Sfr*[A1*v1pr;v1pr'*P*s1pr];
QA=Add(QA,QAfk1,3*(i-1),0);        
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
T1k=-s1pr'*v1pr*Lamk(1)+...
    [(r2'+s2pr'*A2'-r1')*A1*v1pr,v1pr'*A1'*A2*v2pr]*Lamk;
f1k=(T1k/d)-d*Lamk(1);
f2k=-T1k/d;
v12=v1pr'*A1'*(r2d+ph2d*P*A2*s2pr-r1d)/d-ph1d*v1pr'*P*s1pr/d;
[Sfr,Sfrpr]=SfrSfrpr(v12,mus,mud,par);
Sfr=(w/N)*Sfr;
Sfrpr=(w/N)*Sfrpr;
[csf1k,dcsf1k]=csign(f1k,par);
[csf2k,dcsf2k]=csign(f2k,par);
QAfk1=(1/d)*(f1k*csf1k+f2k*csf2k)*Sfr*[A1*v1pr;v1pr'*P*s1pr];
QA=Add(QA,QAfk1,3*(i-1),0);
QAfk2=(1/d)*(f1k*csf1k+f2k*csf2k)*Sfr*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr];
QA=Add(QA,QAfk2,3*(j-1),0);    
end
F(k)=(f1k*csf1k+f2k*csf2k);
end

if PJDT(1,k)==3     %Distance constraint
[i,j,s1pr,s2pr,d,ms,nm]=DistPart(k,PJDT);    
F(k)=Lam(ms);    
end

k=k+1;
end
QAw=QA;

end

