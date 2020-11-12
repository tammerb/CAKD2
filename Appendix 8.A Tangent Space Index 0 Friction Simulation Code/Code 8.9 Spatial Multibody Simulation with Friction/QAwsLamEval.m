function QAwsLam=QAwsLamEval(tn,q,qd,Lam,SJDT,STSDAT,par,w,N)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

QAwsLam=zeros(ngc,nc);

I3=eye(3);
z3=zeros(3,1);
z4=zeros(4,1);

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
Sfr1=(w/N)*Sfr1;
Sfrpr1=(w/N)*Sfrpr1;
F12kcyl=f1prk+f2prk;
f12kcylfr=-F12kcyl*Sfr1;
Ep1=EEval(p1);
Ep2=EEval(p2);
omeg12k=2*vz1pr'*A1'*(Ep1*p1d-Ep2*p2d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
Sfr2=(w/N)*Sfr2;
Sfrpr2=(w/N)*Sfrpr2;
tau12kcylfr=-R*F12kcyl*Sfr2;

%Evaluate QAwsLam
F1prksLam=-A1'*Phisr1k';
T1prksLam=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k');
fx1prksLam=-vx1pr'*F1prksLam+(1/a)*vy1pr'*T1prksLam;
fy1prksLam=-vy1pr'*F1prksLam-(1/a)*vx1pr'*T1prksLam;
fx2prksLam=-(1/a)*vy1pr'*T1prksLam;
fy2prksLam=(1/a)*vx1pr'*T1prksLam;
F12kcylsLam=(1/(f1prk+10^-6))*(fx1prk*fx1prksLam+fy1prk*fy1prksLam)+...
    (1/(f2prk+10^-6))*(fx2prk*fx2prksLam+fy2prk*fy2prksLam);
f12kcylfrsLam=-Sfr1*F12kcylsLam;
tau12kcylfrsLam=-R*Sfr2*F12kcylsLam;
Q1kcylfrsLam=[A1*vz1pr*f12kcylfrsLam;BT1'*A1*vz1pr*f12kcylfrsLam+...
    2*Ep1'*A1*vz1pr*tau12kcylfrsLam];
QAwsLam=Add(QAwsLam,Q1kcylfrsLam,7*(i-1),ms-1);

if j>=1
Q2kcylfrsLam=[-A1*vz1pr*f12kcylfrsLam;-BT2'*A1*vz1pr*f12kcylfrsLam-...
    2*Ep2'*A1*vz1pr*tau12kcylfrsLam];   
QAwsLam=Add(QAwsLam,Q2kcylfrsLam,7*(j-1),ms-1);    
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
Sfr2=(w/N)*Sfr2;
Sfrpr2=(w/N)*Sfrpr2;
tau12kcylfr=-R*F12kcyl*Sfr2;
[csfz1prk,dcsfz1prk]=csign(fz1prk,par);
tau12prkrev5fr=-R*fz1prk*csfz1prk*Sfr2;

%Evaluate QAsLam
F1prksLam=-A1'*Phisr1k';
T1prksLam=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k');
fx1prksLam=-vx1pr'*F1prksLam+(1/a)*vy1pr'*T1prksLam;
fy1prksLam=-vy1pr'*F1prksLam-(1/a)*vx1pr'*T1prksLam;
fx2prksLam=-(1/a)*vy1pr'*T1prksLam;
fy2prksLam=(1/a)*vx1pr'*T1prksLam;
fz1prksLam=vz1pr'*F1prksLam;
F12kcylsLam=(1/(f1prk+10^-6))*(fx1prk*fx1prksLam+fy1prk*fy1prksLam)+...
    (1/(f2prk+10^-6))*(fx2prk*fx2prksLam+fy2prk*fy2prksLam);
tau12kcylfrsLam=-R*Sfr2*F12kcylsLam;
tau12prkrev5frsLam=-R*Sfr2*(csfz1prk+fz1prk*dcsfz1prk)*fz1prksLam;
Q1krevfrsLam=[zeros(3,1);2*Ep1'*A1*vz1pr]*...
    (tau12prkrev5frsLam+tau12kcylfrsLam);
QAwsLam=Add(QAwsLam,Q1krevfrsLam,7*(i-1),ms-1);

if j>=1
Q2krevfrsLam=[zeros(3,1);-2*Ep2'*A1*vz1pr]*...
    (tau12prkrev5frsLam+tau12kcylfrsLam);
QAwsLam=Add(QAwsLam,Q2krevfrsLam,7*(j-1),ms-1);
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
Sfr1=(w/N)*Sfr1;
Sfrpr1=(w/N)*Sfrpr1;
Ep1=EEval(p1);
Ep2=EEval(p2);
f12ktranfr=-F12ktran*Sfr1;


%Evaluate QAsLam
F1prksLam=-A1'*Phisr1k';
T1prksLam=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k');
fx1prksLam=-vx1pr'*F1prksLam+((1/a)*vy1pr'+(1/b)*vz1pr')*T1prksLam;
fy1prksLam=-vy1pr'*F1prksLam-(1/a)*vx1pr'*T1prksLam;
fx2prksLam=-(1/a)*vy1pr'*T1prksLam;
fy2prksLam=(1/a)*vx1pr'*T1prksLam;
fx3prksLam=(1/b)*vz1pr'*T1prksLam;
F12ktransLam=(csfx1prk+fx1prk*dcsfx1prk)*fx1prksLam+...
    (csfy1prk+fy1prk*dcsfy1prk)*fy1prksLam+...
    (csfx2prk+fx2prk*dcsfx2prk)*fx2prksLam+...
    (csfy2prk+fy2prk*dcsfy2prk)*fy2prksLam+...
    (csfx3prk+fx3prk*dcsfx3prk)*fx3prksLam;
f12ktranfrsLam=-Sfr1*F12ktransLam;
Q1ktranfrsLam=[A1*vz1pr;BT1'*A1*vz1pr]*f12ktranfrsLam;
QAwsLam=Add(QAwsLam,Q1ktranfrsLam,7*(i-1),ms-1);

if j>=1
Q2ktranfrsLam=[-A1*vz1pr;-BT2'*A1*vz1pr]*f12ktranfrsLam;  
QAwsLam=Add(QAwsLam,Q2ktranfrsLam,7*(j-1),ms-1);    
end
    
end
k=k+1;

end

end




