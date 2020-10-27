function QA=QAEval(tn,q,qd,SMDT,STSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

uz=[0;0;1];

QA=zeros(ngc,1);

%Account for gravitational force in negative z direction
i=1;
while i<=nb
mi=SMDT(1,i);
QAgi=[-mi*g*uz;zeros(4,1)];
QA=Add(QA,QAgi,7*(i-1),0);
i=i+1;
end

%Account for TSDA forces

TS=1;
while TS<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=STSDATPart(STSDAT,TS);
[r1,p1]=qPart(q,i);
[r1d,p1d]=qPart(qd,i);
r2=[0;0;0];
p2=[1;0;0;0];
r2d=[0;0;0];
p2d=zeros(4,1);
if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
end
A1=ATran(p1);
A2=ATran(p2);
d12=r2+A2*s2pr-r1-A1*s1pr;
BT1=BTran(p1,s1pr);
BT2=BTran(p2,s2pr);
el=sqrt(d12'*d12);
eld=(1/el)*d12'*(r2d+BT2*p2d-r1d-BT1*p1d);
f=K*(el-el0)+C*eld+F;   %User insert F(el,eld) if needed
QA1=(f/el)*[d12;BT1'*d12];
QA=Add(QA,QA1,7*(i-1),0);
if j>=1
QA2=-(f/el)*[d12;BT2'*d12];
QA=Add(QA,QA2,7*(j-1),0);
end

TS=TS+1;
end

end

