function QA=QAEval(tn,q,qd,PMDT,PTSDAT,PRSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA,NRSDA]=...
    parPart(par);

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
el=sqrt(d12'*d12);
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

end

