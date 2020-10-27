function [QAsq,QAsqd]=QAsqqd(tn,q,qd,SMDT,STSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

uy=[0;0;1];

QAsq=zeros(ngc,ngc);
QAsqd=zeros(ngc,ngc);

I3=eye(3);
Z3=zeros(3,1);

T=1;
while T<=NTSDA
    
%Evaluate QAsq
[i,j,s1pr,s2pr,K,C,el0,F]=STSDATPart(STSDAT,T);
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
BT1=BTran(p1,s1pr);
BT2=BTran(p2,s2pr);
BT1d=BTran(p1d,s1pr);
BT2d=BTran(p2d,s2pr);
d12=r2+A2*s2pr-r1-A1*s1pr;
el=sqrt(d12'*d12);
a=r2d+BT2*p2d-r1d-BT1*p1d;
eld=(1/el)*d12'*a;
f=K*(el-el0)+C*eld+F;   %User insert F(el,eld) if needed
elsq1=-(1/el)*d12'*[I3,BT1];
eldsq1=-(1/el^2)*d12'*a*elsq1-(1/el)*a'*[I3,BT1]+...
    (1/el)*d12'*[zeros(3,3),BT1d];
b1=(1/el)*(K*elsq1+C*eldsq1)-(f/el^2)*elsq1;
Q1sq1=[d12;BT1'*d12]*b1+(f/el)*[[-I3,-BT1];...
    [BT1',BT1'*BT1+KEval(s1pr,d12)]];
QAsq=Add(QAsq,Q1sq1,7*(i-1),7*(i-1));

if j>=1
elsq2=(1/el)*d12'*[I3,BT2];
eldsq2=-(1/el^2)*d12'*a*elsq2+(1/el)*a'*[I3,BT2]+...
    (1/el)*d12'*[zeros(3,3),BT2d];
b2=(1/el)*(K*elsq2+C*eldsq2)-(f/el^2)*elsq2;    
Q1sq2=[d12;BT1'*d12]*b2+(f/el)*[[I3,BT2];...
    BT1'*[I3,BT2]];
QAsq=Add(QAsq,Q1sq2,7*(i-1),7*(j-1));
Q2sq1=-[d12;BT2'*d12]*b1-(f/el)*[[-I3,-BT1];...
    BT2'*[-I3,-BT1]];
QAsq=Add(QAsq,Q2sq1,7*(j-1),7*(i-1));
Q2sq2=-[d12;BT2'*d12]*b2-(f/el)*[[I3,BT2];...
    [BT2',BT2'*BT2+KEval(s2pr,d12)]];
QAsq=Add(QAsq,Q2sq2,7*(j-1),7*(j-1));
end

%Evaluate QAsqd
eldsq1d=(1/el)*d12'*[-I3,-BT1];
Q1sq1d=(C/el)*[d12;BT1'*d12]*eldsq1d;
QAsqd=Add(QAsqd,Q1sq1d,7*(i-1),7*(i-1));

if j>=1
eldsq2d=(1/el)*d12'*[I3,BT2];    
Q1sq2d=(C/el)*[d12;BT1'*d12]*eldsq2d;
QAsqd=Add(QAsqd,Q1sq2d,7*(i-1),7*(j-1));
Q2sq1d=-(C/el)*[d12;BT2'*d12]*eldsq1d;
QAsqd=Add(QAsqd,Q2sq1d,7*(j-1),7*(i-1));
Q2sq2d=-(C/el)*[d12;BT2'*d12]*eldsq2d;
QAsqd=Add(QAsqd,Q2sq2d,7*(j-1),7*(j-1));
end   

T=T+1;
end

end


