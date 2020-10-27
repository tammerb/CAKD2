function [QAsq,QAsqd]=QAsqqd(tn,q,qd,Lam,SJDT,STSDAT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

QAsq=zeros(ngc,ngc);
QAsqd=zeros(ngc,ngc);

I3=eye(3);
z3=zeros(3,1);
z4=zeros(4,1);

T=1;
while T<=NTSDA      %Contribution from TSDA
    
%Evaluate QAsq
[i,j,s1pr,s2pr,K,C,el0,F]=STSDATPart(STSDAT,T);
[r1,p1]=qPart(q,i);
[r1d,p1d]=qPart(qd,i);
r2=z3;
p2=[1;0;0;0];
r2d=z3;
p2d=z4;
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

k=1;
while k<=nh    %Contribution from Cylindrical, Revolute, and Translational
                %Joint Friction

if SJDT(1,k)==3     %Cylindrical joint
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm]=CylPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;
Lamk=[Lam(ms);Lam(ms+1);Lam(ms+2);Lam(ms+3)];
[r1,p1]=qPart(q,i);
A1=ATran(p1);
[r1d,p1d]=qPart(qd,i);
BT1=BTran(p1,s1pr);
if j==0
r2=z3;
p2=[1;0;0;0];
r2d=z3;
p2d=z4;    
end
if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
end
BT2=BTran(p2,s2pr);

%q1 Jacobian Calculation
[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41];
Phisr1k=[Phiq1k(:,1),Phiq1k(:,2),Phiq1k(:,3)];
Phisp1k=[Phiq1k(:,4),Phiq1k(:,5),Phiq1k(:,6),Phiq1k(:,7)];

%Reaction Force/Torque Calculaation
F1prk=-A1'*Phisr1k'*Lamk;
T1prk=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k')*Lamk;
fx1prk=-vx1pr'*F1prk+(1/a)*vy1pr'*T1prk;
fy1prk=-vy1pr'*F1prk-(1/a)*vx1pr'*T1prk;
fx2prk=-(1/a)*vy1pr'*T1prk;
fy2prk=(1/a)*vx1pr'*T1prk;
f1prk=sqrt(fx1prk^2+fy1prk^2);
f2prk=sqrt(fx2prk^2+fy2prk^2);

%Friction Force/Torque Calculation
s12dk=-vz1pr'*A1'*(r2d+BT2*p2d-r1d-BT1*p1d);
[Sfr1,Sfrpr1]=SfrSfrpr(s12dk,mus,mud,par);
F12kcyl=f1prk+f2prk;
f12kcylfr=-F12kcyl*Sfr1;
Ep1=EEval(p1);
Ep2=EEval(p2);
omeg12k=2*vz1pr'*A1'*(Ep1*p1d-Ep2*p2d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
tau12kcylfr=-R*F12kcyl*Sfr2;

%Calculate P4 terms in Eq. 6.8.25, Using P4 submatrices
[P4cyl111,P4cyl112,P4cyl122]=...
    bbP4dot2(i,j,vx2pr,s1pr,s2pr,tn,q,Lamk(1),par);
[P4cyl211,P4cyl212,P4cyl222]=...
    bbP4dot2(i,j,vy2pr,s1pr,s2pr,tn,q,Lamk(2),par);
[P4cyl311,P4cyl312,P4cyl322]=bbP4dot1(i,j,vz1pr,vx2pr,tn,q,Lamk(3),par);
[P4cyl411,P4cyl412,P4cyl422]=bbP4dot1(i,j,vz1pr,vy2pr,tn,q,Lamk(4),par);
P411k=P4cyl111+P4cyl211+P4cyl311+P4cyl411;
P412k=P4cyl112+P4cyl212+P4cyl312+P4cyl412;
P41=[P411k,P412k];
Phiksr1TLamsq12=[P41(1,:);P41(2,:);P41(3,:)];
Phiksp1TLamsq12=[P41(4,:);P41(5,:);P41(6,:);P41(7,:)];

%Derivatives of Reaction Force/Torque w.r.t q12
CT1=CTran(p1,Phisr1k'*Lamk);
F1prksq12=-A1'*Phiksr1TLamsq12-[zeros(3,3),CT1,zeros(3,7)];
T1prksq12=-(0.5*GEval(p1)*Phiksp1TLamsq12-atil(s1pr)*A1'*Phiksr1TLamsq12)...
    -[zeros(3,3),-0.5*GEval(Phisp1k'*Lamk)-atil(s1pr)*CT1,zeros(3,7)];
fx1prksq12=-vx1pr'*F1prksq12+(1/a)*vy1pr'*T1prksq12;
fy1prksq12=-vy1pr'*F1prksq12-(1/a)*vx1pr'*T1prksq12;
fx2prksq12=-(1/a)*vy1pr'*T1prksq12;
fy2prksq12=(1/a)*vx1pr'*T1prksq12;
F12kcylsq12=(1/(f1prk+10^-6))*(fx1prk*fx1prksq12+fy1prk*fy1prksq12)+...
    (1/(f2prk+10^-6))*(fx2prk*fx2prksq12+fy2prk*fy2prksq12);

%Derivatives of Friction Force/Torque w.r.t. q12
s12dksq12=-vz1pr'*([zeros(3,3),CTran(p1,(r2d+BT2*p2d-r1d-BT1*p1d)),...
    zeros(3,7)]+A1'*[zeros(3,3),-BTran(p1d,s1pr),zeros(3,3),...
    BTran(p2d,s2pr)]);
omeg12ksq12=2*vz1pr'*[zeros(3,3),CTran(p1,(Ep1*p1d-Ep2*p2d))-...
    A1'*EEval(p1d),zeros(3,3),A1'*EEval(p2d)];
f12kcylfrsq12=-F12kcyl*Sfrpr1*s12dksq12-Sfr1*F12kcylsq12;
tau12kcylfrsq12=-(R^2)*F12kcyl*Sfrpr2*omeg12ksq12-R*Sfr2*F12kcylsq12;

%Evaluate QAsq12
a1=[zeros(3,3),BTran(p1,vz1pr*f12kcylfr),zeros(3,7)];
a2=[zeros(4,3),KEval(s1pr,A1*vz1pr*f12kcylfr)+...
    BT1'*BTran(p1,vz1pr*f12kcylfr),zeros(4,7)]+...
    2*[zeros(4,3),REval(A1*vz1pr*tau12kcylfr)+...
    Ep1'*BTran(p1,vz1pr*tau12kcylfr),zeros(4,7)];
Q1kcylfrsq12=[A1*vz1pr*f12kcylfrsq12+a1;BT1'*A1*vz1pr*f12kcylfrsq12+...
    2*Ep1'*A1*vz1pr*tau12kcylfrsq12+a2];
Q1kcylfrsq1=[Q1kcylfrsq12(:,1),Q1kcylfrsq12(:,2),Q1kcylfrsq12(:,3),...
    Q1kcylfrsq12(:,4),Q1kcylfrsq12(:,5),Q1kcylfrsq12(:,6),...
    Q1kcylfrsq12(:,7)];
QAsq=Add(QAsq,Q1kcylfrsq1,7*(i-1),7*(i-1));

if j>=1
Q1kcylfrsq2=[Q1kcylfrsq12(:,8),Q1kcylfrsq12(:,9),Q1kcylfrsq12(:,10),...
    Q1kcylfrsq12(:,11),Q1kcylfrsq12(:,12),Q1kcylfrsq12(:,13),...
    Q1kcylfrsq12(:,14)];
QAsq=Add(QAsq,Q1kcylfrsq2,7*(i-1),7*(j-1));
a3=[zeros(4,3),BT2'*BTran(p1,vz1pr*f12kcylfr),zeros(4,3),...
    KEval(s2pr,A1*vz1pr*f12kcylfr)]+...
    2*[zeros(4,3),Ep2'*BTran(p1,vz1pr*tau12kcylfr),zeros(4,3),...
    REval(A1*vz1pr*tau12kcylfr)];
Q2kcylfrsq12=[-A1*vz1pr*f12kcylfrsq12-a1;-BT2'*A1*vz1pr*f12kcylfrsq12-...
    2*Ep2'*A1*vz1pr*tau12kcylfrsq12-a3];
Q2kcylfrsq1=[Q2kcylfrsq12(:,1),Q2kcylfrsq12(:,2),Q2kcylfrsq12(:,3),...
    Q2kcylfrsq12(:,4),Q2kcylfrsq12(:,5),Q2kcylfrsq12(:,6),...
    Q2kcylfrsq12(:,7)];
Q2kcylfrsq2=[Q2kcylfrsq12(:,8),Q2kcylfrsq12(:,9),Q2kcylfrsq12(:,10),...
    Q2kcylfrsq12(:,11),Q2kcylfrsq12(:,12),Q2kcylfrsq12(:,13),...
    Q2kcylfrsq12(:,14)];
QAsq=Add(QAsq,Q2kcylfrsq1,7*(j-1),7*(i-1));
QAsq=Add(QAsq,Q2kcylfrsq2,7*(j-1),7*(j-1));   
end

%Evaluate QAsq12d
s12dksq12d=-vz1pr'*A1'*[-I3,-BT1,I3,BT2];
omeg12ksq12d=2*vz1pr'*A1'*[zeros(3,3),Ep1,zeros(3,3),-Ep2];
f12kcylfrsq12d=-F12kcyl*Sfrpr1*s12dksq12d;
tau12kcylfrsq12d=-(R^2)*F12kcyl*Sfrpr2*omeg12ksq12d;
Q1kcylfrsq12d=[A1*vz1pr*f12kcylfrsq12d;BT1'*A1*vz1pr*f12kcylfrsq12d+...
    2*Ep1'*A1*vz1pr*tau12kcylfrsq12d];
Q1kcylfrsq1d=[Q1kcylfrsq12d(:,1),Q1kcylfrsq12d(:,2),Q1kcylfrsq12d(:,3),...
    Q1kcylfrsq12d(:,4),Q1kcylfrsq12d(:,5),Q1kcylfrsq12d(:,6),...
    Q1kcylfrsq12d(:,7)];
QAsqd=Add(QAsqd,Q1kcylfrsq1d,7*(i-1),7*(i-1));

if j>=1
Q1kcylfrsq2d=[Q1kcylfrsq12d(:,8),Q1kcylfrsq12d(:,9),Q1kcylfrsq12d(:,10),...
    Q1kcylfrsq12d(:,11),Q1kcylfrsq12d(:,12),Q1kcylfrsq12d(:,13),...
    Q1kcylfrsq12d(:,14)];
QAsqd=Add(QAsqd,Q1kcylfrsq2d,7*(i-1),7*(j-1));
Q2kcylfrsq12d=[-A1*vz1pr*f12kcylfrsq12d;-BT2'*A1*vz1pr*f12kcylfrsq12d-...
    2*Ep2'*A1*vz1pr*tau12kcylfrsq12d];
Q2kcylfrsq1d=[Q2kcylfrsq12d(:,1),Q2kcylfrsq12d(:,2),Q2kcylfrsq12d(:,3),...
    Q2kcylfrsq12d(:,4),Q2kcylfrsq12d(:,5),Q2kcylfrsq12d(:,6),...
    Q2kcylfrsq12d(:,7)];
Q2kcylfrsq2d=[Q2kcylfrsq12d(:,8),Q2kcylfrsq12d(:,9),Q2kcylfrsq12d(:,10),...
    Q2kcylfrsq12d(:,11),Q2kcylfrsq12d(:,12),Q2kcylfrsq12d(:,13),...
    Q2kcylfrsq12d(:,14)];
QAsqd=Add(QAsqd,Q2kcylfrsq1d,7*(j-1),7*(i-1));
QAsqd=Add(QAsqd,Q2kcylfrsq2d,7*(j-1),7*(j-1));
end

end     %End Cylindrical Joint Contributions

if SJDT(1,k)==4     %Revolute joint
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm]=RevPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;
Lamk=[Lam(ms);Lam(ms+1);Lam(ms+2);Lam(ms+3);Lam(ms+4)];
[r1,p1]=qPart(q,i);
A1=ATran(p1);
[r1d,p1d]=qPart(qd,i);
if j==0
r2=z3;
p2=[1;0;0;0];
r2d=z3;
p2d=z4;
end
if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
end

%q1 Jacobian Calcualtion
[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
[Phiqrev51,Phiqrev52]=bbPhiqdot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);
Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41;Phiqrev51];
Phisr1k=[Phiq1k(:,1),Phiq1k(:,2),Phiq1k(:,3)];
Phisp1k=[Phiq1k(:,4),Phiq1k(:,5),Phiq1k(:,6),Phiq1k(:,7)];

%Reaction Force/Torque Calculation
F1prk=-A1'*Phisr1k'*Lamk;
T1prk=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k')*Lamk;
fx1prk=-vx1pr'*F1prk+(1/a)*vy1pr'*T1prk;
fy1prk=-vy1pr'*F1prk-(1/a)*vx1pr'*T1prk;
fx2prk=-(1/a)*vy1pr'*T1prk;
fy2prk=(1/a)*vx1pr'*T1prk;
fz1prk=vz1pr'*F1prk;
f1prk=sqrt(fx1prk^2+fy1prk^2);
f2prk=sqrt(fx2prk^2+fy2prk^2);
F12kcyl=f1prk+f2prk;

%Friction Force/Torque Calculation
Ep1=EEval(p1);
Ep2=EEval(p2);
BT1=BTran(p1,s1pr);
BT2=BTran(p2,s2pr);
BTv1=BTran(p1,vz1pr);
omeg12k=2*vz1pr'*A1'*(Ep1*p1d-Ep2*p2d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
tau12kcylfr=-R*F12kcyl*Sfr2;
[csfz1prk,dcsfz1prk]=csign(fz1prk,par);
tau12prkrev5fr=-R*fz1prk*csfz1prk*Sfr2;

%Calculate P4 terms in Eq. 6.8.25, Using P4 submatrices
[P4cyl111,P4cyl112,P4cyl122]=...
    bbP4dot2(i,j,vx2pr,s1pr,s2pr,tn,q,Lamk(1),par);
[P4cyl211,P4cyl212,P4cyl222]=...
    bbP4dot2(i,j,vy2pr,s1pr,s2pr,tn,q,Lamk(2),par);
[P4cyl311,P4cyl312,P4cyl322]=bbP4dot1(i,j,vz1pr,vx2pr,tn,q,Lamk(3),par);
[P4cyl411,P4cyl412,P4cyl422]=bbP4dot1(i,j,vz1pr,vy2pr,tn,q,Lamk(4),par);
[P4rev511,P4rev512,P4rev522]=...
    bbP4dot2(i,j,vz2pr,s1pr,s2pr,tn,q,Lamk(4),par);
P411=P4cyl111+P4cyl211+P4cyl311+P4cyl411+P4rev511;
P412=P4cyl112+P4cyl212+P4cyl312+P4cyl412+P4rev512;
P41=[P411,P412];
Phiksr1TLamsq12=[P41(1,:);P41(2,:);P41(3,:)];
Phiksp1TLamsq12=[P41(4,:);P41(5,:);P41(6,:);P41(7,:)];

%Derivatives of friction Torque w.r.t q12
CT1=CTran(p1,Phisr1k'*Lamk);
F1prksq12=-A1'*Phiksr1TLamsq12-[zeros(3,3),CT1,zeros(3,7)];
T1prksq12=-(0.5*GEval(p1)*Phiksp1TLamsq12-atil(s1pr)*A1'*Phiksr1TLamsq12)...
    -[zeros(3,3),-0.5*GEval(Phisp1k'*Lamk)-atil(s1pr)*CT1,zeros(3,7)];
fx1prksq12=-vx1pr'*F1prksq12+(1/a)*vy1pr'*T1prksq12;
fy1prksq12=-vy1pr'*F1prksq12-(1/a)*vx1pr'*T1prksq12;
fx2prksq12=-(1/a)*vy1pr'*T1prksq12;
fy2prksq12=(1/a)*vx1pr'*T1prksq12;
fz1prksq12=vz1pr'*F1prksq12;
F12kcylsq12=(1/(f1prk+10^-6))*(fx1prk*fx1prksq12+fy1prk*fy1prksq12)+...
    (1/(f2prk+10^-6))*(fx2prk*fx2prksq12+fy2prk*fy2prksq12);
omeg12ksq12=2*vz1pr'*[zeros(3,3),CTran(p1,(Ep1*p1d-Ep2*p2d))-...
    A1'*EEval(p1d),zeros(3,3),A1'*EEval(p2d)];
tau12kcylfrsq12=-(R^2)*F12kcyl*Sfrpr2*omeg12ksq12-R*Sfr2*F12kcylsq12;
tau12prkrev5frsq12=-(R^2)*fz1prk*csfz1prk*Sfrpr2*omeg12ksq12-...
    R*Sfr2*(csfz1prk+fz1prk*dcsfz1prk)*fz1prksq12;

%Evaluate QAsq12
b1=2*[zeros(4,3),REval(A1*vz1pr)+Ep1'*BTv1,zeros(4,7)]*...
    (tau12kcylfr+tau12prkrev5fr);
Q1krevfrsq12=[zeros(3,14);...
    2*Ep1'*A1*vz1pr*(tau12kcylfrsq12+tau12prkrev5frsq12)+b1];
Q1krevfrsq1=[Q1krevfrsq12(:,1),Q1krevfrsq12(:,2),Q1krevfrsq12(:,3),...
    Q1krevfrsq12(:,4),Q1krevfrsq12(:,5),Q1krevfrsq12(:,6),...
    Q1krevfrsq12(:,7)];
QAsq=Add(QAsq,Q1krevfrsq1,7*(i-1),7*(i-1));

if j>=1
Q1krevfrsq2=[Q1krevfrsq12(:,8),Q1krevfrsq12(:,9),Q1krevfrsq12(:,10),...
    Q1krevfrsq12(:,11),Q1krevfrsq12(:,12),Q1krevfrsq12(:,13),...
    Q1krevfrsq12(:,14)];
QAsq=Add(QAsq,Q1krevfrsq2,7*(i-1),7*(j-1)); 
b2=2*[zeros(4,3),Ep2'*BTv1,zeros(4,3),REval(A1*vz1pr)]*...
    (tau12kcylfr+tau12prkrev5fr);
Q2krevfrsq12=[zeros(3,14);...
    -2*Ep2'*A1*vz1pr*(tau12kcylfrsq12+tau12prkrev5frsq12)-b2];
Q2krevfrsq1=[Q2krevfrsq12(:,1),Q2krevfrsq12(:,2),Q2krevfrsq12(:,3),...
    Q2krevfrsq12(:,4),Q2krevfrsq12(:,5),Q2krevfrsq12(:,6),...
    Q2krevfrsq12(:,7)];
Q2krevfrsq2=[Q2krevfrsq12(:,8),Q2krevfrsq12(:,9),Q2krevfrsq12(:,10),...
    Q2krevfrsq12(:,11),Q2krevfrsq12(:,12),Q2krevfrsq12(:,13),...
    Q2krevfrsq12(:,14)];
QAsq=Add(QAsq,Q2krevfrsq1,7*(j-1),7*(i-1));
QAsq=Add(QAsq,Q2krevfrsq2,7*(j-1),7*(j-1));    
end

%Evaluate QAsq12d
omeg12ksq12d=2*vz1pr'*A1'*[zeros(3,3),Ep1,zeros(3,3),-Ep2];
tau12kcylfrsq12d=-(R^2)*F12kcyl*Sfrpr2*omeg12ksq12d;
tau12prkrev5frsq12d=-(R^2)*fz1prk*csfz1prk*Sfrpr2*omeg12ksq12d;
Q1krevfrsq12d=[zeros(3,14);...
    2*Ep1'*A1*vz1pr*(tau12kcylfrsq12d+tau12prkrev5frsq12d)];
Q1krevfrsq1d=[Q1krevfrsq12d(:,1),Q1krevfrsq12d(:,2),Q1krevfrsq12d(:,3),...
    Q1krevfrsq12d(:,4),Q1krevfrsq12d(:,5),Q1krevfrsq12d(:,6),...
    Q1krevfrsq12d(:,7)];
QAsqd=Add(QAsqd,Q1krevfrsq1d,7*(i-1),7*(i-1));

if j>=1
Q1krevfrsq2d=[Q1krevfrsq12d(:,8),Q1krevfrsq12d(:,9),Q1krevfrsq12d(:,10),...
    Q1krevfrsq12d(:,11),Q1krevfrsq12d(:,12),Q1krevfrsq12d(:,13),...
    Q1krevfrsq12d(:,14)];
QAsqd=Add(QAsqd,Q1krevfrsq2d,7*(i-1),7*(j-1));
Q2krevfrsq12d=[zeros(3,14);...
    -2*Ep2'*A1*vz1pr*(tau12kcylfrsq12d+tau12prkrev5frsq12d)];
Q2krevfrsq1d=[Q2krevfrsq12d(:,1),Q2krevfrsq12d(:,2),Q2krevfrsq12d(:,3),...
    Q2krevfrsq12d(:,4),Q2krevfrsq12d(:,5),Q2krevfrsq12d(:,6),...
    Q2krevfrsq12d(:,7)];
Q2krevfrsq2d=[Q2krevfrsq12d(:,8),Q2krevfrsq12d(:,9),Q2krevfrsq12d(:,10),...
    Q2krevfrsq12d(:,11),Q2krevfrsq12d(:,12),Q2krevfrsq12d(:,13),...
    Q2krevfrsq12d(:,14)];
QAsqd=Add(QAsqd,Q2krevfrsq1d,7*(j-1),7*(i-1));
QAsqd=Add(QAsqd,Q2krevfrsq2d,7*(j-1),7*(j-1));
end

end

if SJDT(1,k)==5     %Translational joint
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,b,mus,mud,ms,nm]=TranPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;
Lamk=[Lam(ms);Lam(ms+1);Lam(ms+2);Lam(ms+3);Lam(ms+4)];
[r1,p1]=qPart(q,i);
A1=ATran(p1);
[r1d,p1d]=qPart(qd,i);
if j==0
r2=z3;
p2=[1;0;0;0];
r2d=z3;
p2d=z4;
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
end
%q1 Jacobian Calculation
[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
[Phiqtran51,Phiqtran52]=bbPhiqdot1(i,j,vy1pr,vx2pr,tn,q,par);
Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41;Phiqtran51];
Phisr1k=[Phiq1k(:,1),Phiq1k(:,2),Phiq1k(:,3)];
Phisp1k=[Phiq1k(:,4),Phiq1k(:,5),Phiq1k(:,6),Phiq1k(:,7)];

%Reaction Force/Torque Calculation
F1prk=-A1'*Phisr1k'*Lamk;
T1prk=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k')*Lamk;
fx1prk=-vx1pr'*F1prk+((1/a)*vy1pr'+(1/b)*vz1pr')*T1prk;
fy1prk=-vy1pr'*F1prk-(1/a)*vx1pr'*T1prk;
fx2prk=-(1/a)*vy1pr'*T1prk;
fy2prk=(1/a)*vx1pr'*T1prk;
fx3prk=(1/b)*vz1pr'*T1prk;
[csfx1prk,dcsfx1prk]=csign(fx1prk,par);
[csfy1prk,dcsfy1prk]=csign(fy1prk,par);
[csfx2prk,dcsfx2prk]=csign(fx2prk,par);
[csfy2prk,dcsfy2prk]=csign(a,par);
[csfx3prk,dcsfx3prk]=csign(fx3prk,par);
F12ktran=fx1prk*csfx1prk+fy1prk*csfy1prk+fx2prk*csfx2prk+...
    fy2prk*csfy2prk+fx3prk*csfx3prk;

%Friction Force Calculation
BT1=BTran(p1,s1pr);
BT2=BTran(p2,s2pr);
s12dk=-vz1pr'*A1'*(r2d+BT2*p2d-r1d-BT1*p1d);
[Sfr1,Sfrpr1]=SfrSfrpr(s12dk,mus,mud,par);
Ep1=EEval(p1);
Ep2=EEval(p2);
f12ktranfr=-F12ktran*Sfr1;

%Calculate P4 terms in Eq. 6.8.25, Using P4 submatrices
[P4cyl111,P4cyl112,P4cyl122]=...
    bbP4dot2(i,j,vx2pr,s1pr,s2pr,tn,q,Lamk(1),par);
[P4cyl211,P4cyl212,P4cyl222]=...
    bbP4dot2(i,j,vy2pr,s1pr,s2pr,tn,q,Lamk(2),par);
[P4cyl311,P4cyl312,P4cyl322]=bbP4dot1(i,j,vz1pr,vx2pr,tn,q,Lamk(3),par);
[P4cyl411,P4cyl412,P4cyl422]=bbP4dot1(i,j,vz1pr,vy2pr,tn,q,Lamk(4),par);
[P4tran511,P4tran512,P4tran522]=...
    bbP4dot1(i,j,vy1pr,vx2pr,tn,q,Lamk(5),par);
P411k=P4cyl111+P4cyl211+P4cyl311+P4cyl411+P4tran511;
P412k=P4cyl112+P4cyl212+P4cyl312+P4cyl412+P4tran512;
P41=[P411k,P412k];
Phiksr1TLamsq12=[P41(1,:);P41(2,:);P41(3,:)];
Phiksp1TLamsq12=[P41(4,:);P41(5,:);P41(6,:);P41(7,:)];

%Derivatives of Reaction Force w.r.t. q12
CT1=CTran(p1,Phisr1k'*Lamk);
F1prksq12=-A1'*Phiksr1TLamsq12-[zeros(3,3),CT1,zeros(3,7)];
T1prksq12=-(0.5*GEval(p1)*Phiksp1TLamsq12-atil(s1pr)*A1'*Phiksr1TLamsq12)...
    -[zeros(3,3),-0.5*GEval(Phisp1k'*Lamk)-atil(s1pr)*CT1,zeros(3,7)];
fx1prksq12=-vx1pr'*F1prksq12+((1/a)*vy1pr'+(1/b)*vz1pr')*T1prksq12;
fy1prksq12=-vy1pr'*F1prksq12-(1/a)*vx1pr'*T1prksq12;
fx2prksq12=-(1/a)*vy1pr'*T1prksq12;
fy2prksq12=(1/a)*vx1pr'*T1prksq12;
fx3prksq12=(1/b)*vz1pr'*T1prksq12;
F12ktransq12=(csfx1prk+fx1prk*dcsfx1prk)*fx1prksq12+...
    (csfy1prk+fy1prk*dcsfy1prk)*fy1prksq12+...
    (csfx2prk+fx2prk*dcsfx2prk)*fx2prksq12+...
    (csfy2prk+fy2prk*dcsfy2prk)*fy2prksq12+...
    (csfx3prk+fx3prk*dcsfx3prk)*fx3prksq12;

%Evaluate QAsq12
s12dksq12=-vz1pr'*([zeros(3,3),CTran(p1,(r2d+BT2*p2d-r1d-BT1*p1d)),...
    zeros(3,7)]+A1'*[zeros(3,3),-BTran(p1d,s1pr),zeros(3,3),...
    BTran(p2d,s2pr)]);
f12ktranfrsq12=-F12ktran*Sfrpr1*s12dksq12-Sfr1*F12ktransq12;
BT1vz=BTran(p1,vz1pr);
c1=[zeros(3,3),BT1vz,zeros(3,7)]*f12ktranfr;
c2=[zeros(4,3),KEval(s1pr,A1*vz1pr)+BT1'*BT1vz,zeros(4,7)]*f12ktranfr;
Q1ktranfrsq12=[A1*vz1pr*f12ktranfrsq12+c1;BT1'*A1*vz1pr*f12ktranfrsq12+c2];
Q1ktranfrsq1=[Q1ktranfrsq12(:,1),Q1ktranfrsq12(:,2),Q1ktranfrsq12(:,3),...
    Q1ktranfrsq12(:,4),Q1ktranfrsq12(:,5),Q1ktranfrsq12(:,6),...
    Q1ktranfrsq12(:,7)];
QAsq=Add(QAsq,Q1ktranfrsq1,7*(i-1),7*(i-1));

if j>=1
Q1ktranfrsq2=[Q1ktranfrsq12(:,8),Q1ktranfrsq12(:,9),Q1ktranfrsq12(:,10),...
    Q1ktranfrsq12(:,11),Q1ktranfrsq12(:,12),Q1ktranfrsq12(:,13),...
    Q1ktranfrsq12(:,14)];
QAsq=Add(QAsq,Q1ktranfrsq2,7*(i-1),7*(j-1));
c3=[zeros(4,3),BT2'*BT1vz,zeros(4,3),KEval(s2pr,A1*vz1pr)]*f12ktranfr;
Q2ktranfrsq12=[-A1*vz1pr*f12ktranfrsq12-c1;...
    -BT2'*A1*vz1pr*f12ktranfrsq12-c3];
Q2ktranfrsq1=[Q2ktranfrsq12(:,1),Q2ktranfrsq12(:,2),Q2ktranfrsq12(:,3),...
    Q2ktranfrsq12(:,4),Q2ktranfrsq12(:,5),Q2ktranfrsq12(:,6),...
    Q2ktranfrsq12(:,7)];
Q2ktranfrsq2=[Q2ktranfrsq12(:,8),Q2ktranfrsq12(:,9),Q2ktranfrsq12(:,10),...
    Q2ktranfrsq12(:,11),Q2ktranfrsq12(:,12),Q2ktranfrsq12(:,13),...
    Q2ktranfrsq12(:,14)];
QAsq=Add(QAsq,Q2ktranfrsq1,7*(j-1),7*(i-1));
QAsq=Add(QAsq,Q2ktranfrsq2,7*(j-1),7*(j-1));   
end

%Evaluate QAsq12d
s12dksq12d=-vz1pr'*A1'*[-I3,-BT1,I3,BT2];
f12ktranfrsq12d=-F12ktran*Sfrpr1*s12dksq12d;
Q1ktranfrsq12d=[A1*vz1pr;BT1'*A1*vz1pr]*f12ktranfrsq12d;
Q1ktranfrsq1d=[Q1ktranfrsq12d(:,1),Q1ktranfrsq12d(:,2),Q1ktranfrsq12d(:,3),...
    Q1ktranfrsq12d(:,4),Q1ktranfrsq12d(:,5),Q1ktranfrsq12d(:,6),...
    Q1ktranfrsq12d(:,7)];
QAsqd=Add(QAsqd,Q1ktranfrsq1d,7*(i-1),7*(i-1));

if j>=1
Q1ktranfrsq2d=[Q1ktranfrsq12d(:,8),Q1ktranfrsq12d(:,9),...
    Q1ktranfrsq12d(:,10),Q1ktranfrsq12d(:,11),Q1ktranfrsq12d(:,12),...
    Q1ktranfrsq12d(:,13),Q1ktranfrsq12d(:,14)];
QAsqd=Add(QAsqd,Q1ktranfrsq2d,7*(i-1),7*(j-1));
Q2ktranfrsq12d=[-A1*vz1pr;-BT2'*A1*vz1pr]*f12ktranfrsq12d;
Q2ktranfrsq1d=[Q2ktranfrsq12d(:,1),Q2ktranfrsq12d(:,2),...
    Q2ktranfrsq12d(:,3),Q2ktranfrsq12d(:,4),Q2ktranfrsq12d(:,5),...
    Q2ktranfrsq12d(:,6),Q2ktranfrsq12d(:,7)];
Q2ktranfrsq2d=[Q2ktranfrsq12d(:,8),Q2ktranfrsq12d(:,9),...
    Q2ktranfrsq12d(:,10),Q2ktranfrsq12d(:,11),Q2ktranfrsq12d(:,12),...
    Q2ktranfrsq12d(:,13),Q2ktranfrsq12d(:,14)];
QAsqd=Add(QAsqd,Q2ktranfrsq1d,7*(j-1),7*(i-1));
QAsqd=Add(QAsqd,Q2ktranfrsq2d,7*(j-1),7*(j-1));
end
    
end
k=k+1;

end

end



