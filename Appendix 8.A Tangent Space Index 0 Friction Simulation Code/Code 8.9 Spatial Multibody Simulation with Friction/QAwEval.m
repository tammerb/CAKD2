function QA=QAwEval(tn,q,qd,Lam,SMDT,SJDT,STSDAT,par,w,N)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

uz=[0;0;1];
z3=[0;0;0];
z4=[0;0;0;0];
I3=eye(3);

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

T=1;
while T<=NTSDA
[i,j,s1pr,s2pr,K,C,el0,F]=STSDATPart(STSDAT,T);
[r1,p1]=qPart(q,i);
[r1d,p1d]=qPart(qd,i);
r2=[0;0;0];
p2=[1;0;0;0];
r2d=[0;0;0];
p2d=z4;
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

T=T+1;
end

%Account for friction in cylindrical, revolute, and translational joints 
k=1;
while k<=nh
    
if SJDT(1,k)==3     %Cylindrical joint
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm]=CylPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;
[r1,p1]=qPart(q,i);
A1=ATran(p1);
[r1d,p1d]=qPart(qd,i);
Ep1=EEval(p1);
BT1=BTran(p1,s1pr);
Lamk=[Lam(ms);Lam(ms+1);Lam(ms+2);Lam(ms+3)];
[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41];
Phisr1k=[Phiq1k(:,1),Phiq1k(:,2),Phiq1k(:,3)];
Phisp1k=[Phiq1k(:,4),Phiq1k(:,5),Phiq1k(:,6),Phiq1k(:,7)];
F1prk=-A1'*Phisr1k'*Lamk;
T1prk=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k')*Lamk;
fx1prk=-vx1pr'*F1prk+(1/a)*vy1pr'*T1prk;
fy1prk=-vy1pr'*F1prk-(1/a)*vx1pr'*T1prk;
fx2prk=-(1/a)*vy1pr'*T1prk;
fy2prk=(1/a)*vx1pr'*T1prk;
f1prk=sqrt(fx1prk^2+fy1prk^2);
f2prk=sqrt(fx2prk^2+fy2prk^2);
F12kcyl=f1prk+f2prk;

if j==0
s12dk=-vz1pr'*A1'*(-r1d-BT1*p1d);
[Sfr1,Sfrpr1]=SfrSfrpr(s12dk,mus,mud,par);
Sfr1=(w/N)*Sfr1;
Sfrpr1=(w/N)*Sfrpr1;
f12kcylfr=-F12kcyl*Sfr1;
omeg12k=2*vz1pr'*A1'*(Ep1*p1d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
Sfr2=(w/N)*Sfr2;
Sfrpr2=(w/N)*Sfrpr1;
tau12kcylfr=-R*F12kcyl*Sfr2;
Q1kcylfr=[A1*vz1pr*f12kcylfr;BT1'*A1*vz1pr*f12kcylfr+...
    2*Ep1'*A1*vz1pr*tau12kcylfr];
QA=Add(QA,Q1kcylfr,7*(i-1),0);    
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
Ep2=EEval(p2);
BT2=BTran(p2,s2pr);
s12dk=-vz1pr'*A1'*(r2d+BT2*p2d-r1d-BT1*p1d);
[Sfr1,Sfrpr1]=SfrSfrpr(s12dk,mus,mud,par);
Sfr1=(w/N)*Sfr1;
Sfrpr1=(w/N)*Sfrpr1;
f12kcylfr=-F12kcyl*Sfr1;
omeg12k=2*vz1pr'*A1'*(Ep1*p1d-Ep2*p2d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
Sfr2=(w/N)*Sfr2;
Sfrpr2=(w/N)*Sfrpr2;
tau12kcylfr=-R*F12kcyl*Sfr2;
Q1kcylfr=[A1*vz1pr*f12kcylfr;BT1'*A1*vz1pr*f12kcylfr+...
    2*Ep1'*A1*vz1pr*tau12kcylfr];
QA=Add(QA,Q1kcylfr,7*(i-1),0);
Q2kcylfr=[-A1*vz1pr*f12kcylfr;-BT2'*A1*vz1pr*f12kcylfr-...
    2*Ep2'*A1*vz1pr*tau12kcylfr];
QA=Add(QA,Q2kcylfr,7*(j-1),0);    
end
end

if SJDT(1,k)==4     %Revolute joint
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm]=RevPart(k,SJDT);
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;
Lamk=[Lam(ms);Lam(ms+1);Lam(ms+2);Lam(ms+3);Lam(ms+4)];
[r1,p1]=qPart(q,i);
A1=ATran(p1);
[r1d,p1d]=qPart(qd,i);
Ep1=EEval(p1);
[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
[Phiqrev51,Phiqrev52]=bbPhiqdot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);
Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41;Phiqrev51];
Phisr1k=[Phiq1k(:,1),Phiq1k(:,2),Phiq1k(:,3)];
Phisp1k=[Phiq1k(:,4),Phiq1k(:,5),Phiq1k(:,6),Phiq1k(:,7)];
F1prk=-A1'*Phisr1k'*Lamk;
T1prk=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k')*Lamk;
fx1prk=-vx1pr'*F1prk+(1/a)*vy1pr'*T1prk;
fy1prk=-vy1pr'*F1prk-(1/a)*vx1pr'*T1prk;
fx2prk=-(1/a)*vy1pr'*T1prk;
fy2prk=(1/a)*vx1pr'*T1prk;
f1prk=sqrt(fx1prk^2+fy1prk^2);
f2prk=sqrt(fx2prk^2+fy2prk^2);
F12kcyl=f1prk+f2prk;

if j==0
omeg12k=2*vz1pr'*A1'*(Ep1*p1d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
Sfr2=(w/N)*Sfr2;
Sfrpr2=(w/N)*Sfrpr2;
tau12kcylfr=-R*F12kcyl*Sfr2;
fz1prk=vz1pr'*F1prk;
[csfz1prk,dcsfz1prk]=csign(fz1prk,par);
tau12prkrev5fr=-R*fz1prk*csfz1prk*Sfr2;
Q1krevfr=[0;0;0;2*Ep1'*A1*vz1pr*(tau12kcylfr+tau12prkrev5fr)];
QA=Add(QA,Q1krevfr,7*(i-1),0);    
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
Ep2=EEval(p2);
omeg12k=2*vz1pr'*A1'*(Ep1*p1d-Ep2*p2d);
[Sfr2,Sfrpr2]=SfrSfrpr(R*omeg12k,mus,mud,par);
Sfr2=(w/N)*Sfr2;
Sfrpr2=(w/N)*Sfrpr2;
tau12kcylfr=-R*F12kcyl*Sfr2;
fz1prk=vz1pr'*F1prk;
[csfz1prk,dcsfz1prk]=csign(fz1prk,par);
tau12prkrev5fr=-R*fz1prk*csfz1prk*Sfr2;
Q1krevfr=[0;0;0;2*Ep1'*A1*vz1pr*(tau12kcylfr+tau12prkrev5fr)];
QA=Add(QA,Q1krevfr,7*(i-1),0);
Q2krevfr=[0;0;0;-2*Ep2'*A1*vz1pr*(tau12kcylfr+tau12prkrev5fr)];
QA=Add(QA,Q2krevfr,7*(j-1),0);    
end
end

if SJDT(1,k)==5     %Translational joint
[i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,b,mus,mud,ms,nm]=TranPart(k,SJDT);    
vy1pr=atil(vz1pr)*vx1pr;   
vy2pr=atil(vz2pr)*vx2pr;
[r1,p1]=qPart(q,i);
A1=ATran(p1);
[r1d,p1d]=qPart(qd,i);
Lamk=[Lam(ms);Lam(ms+1);Lam(ms+2);Lam(ms+3);Lam(ms+4)];
[Phiqcyl11,Phiqcyl12]=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl21,Phiqcyl22]=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
[Phiqcyl31,Phiqcyl32]=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
[Phiqcyl41,Phiqcyl42]=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
[Phiqtran51,Phiqtran52]=bbPhiqdot1(i,j,vy1pr,vx2pr,tn,q,par);
Phiq1k=[Phiqcyl11;Phiqcyl21;Phiqcyl31;Phiqcyl41;Phiqtran51];
Phisr1k=[Phiq1k(:,1),Phiq1k(:,2),Phiq1k(:,3)];
Phisp1k=[Phiq1k(:,4),Phiq1k(:,5),Phiq1k(:,6),Phiq1k(:,7)];
F1prk=-A1'*Phisr1k'*Lamk;
T1prk=-(0.5*GEval(p1)*Phisp1k'-atil(s1pr)*A1'*Phisr1k')*Lamk;
fx1prk=-vx1pr'*F1prk+((1/a)*vy1pr'+(1/b)*vz1pr')*T1prk;
fy1prk=-vy1pr'*F1prk-(1/a)*vx1pr'*T1prk;
fx2prk=-(1/a)*vy1pr'*T1prk;
fy2prk=(1/a)*vx1pr'*T1prk;
f3prk=(1/b)*vz1pr'*T1prk;
[csfx1prk,dcsfx1prk]=csign(fx1prk,par);
[csfy1prk,dcsfy1prk]=csign(fy1prk,par);
[csfx2prk,dcsfx2prk]=csign(fx2prk,par);
[csfy2prk,dcsfy2prk]=csign(fy2prk,par);
[csf3prk,dcsf3prk]=csign(f3prk,par);
F12ktran=fx1prk*csfx1prk+fy1prk*csfy1prk+fx2prk*csfx2prk+...
    fy2prk*csfy2prk+f3prk*csf3prk;
BT1=BTran(p1,s1pr);

if j==0
s12dk=-vz1pr'*A1'*(-r1d-BT1*p1d);
[Sfr1,Sfrpr1]=SfrSfrpr(s12dk,mus,mud,par);
Sfr1=(w/N)*Sfr1;
Sfrpr1=(w/N)*Sfrpr1;
f12ktranfr=-F12ktran*Sfr1;
Q1ktranfr=[A1*vz1pr;BT1'*A1*vz1pr]*f12ktranfr;
QA=Add(QA,Q1ktranfr,7*(i-1),0);        
end

if j>=1
[r2,p2]=qPart(q,j);
[r2d,p2d]=qPart(qd,j);
BT2=BTran(p2,s2pr);
s12dk=-vz1pr'*A1'*(r2d+BT2*p2d-r1d-BT1*p1d);
[Sfr1,Sfrpr1]=SfrSfrpr(s12dk,mus,mud,par);
Sfr1=(w/N)*Sfr1;
Sfrpr1=(w/N)*Sfrpr1;
f12ktranfr=-F12ktran*Sfr1;
Q1ktranfr=[A1*vz1pr;BT1'*A1*vz1pr]*f12ktranfr;
QA=Add(QA,Q1ktranfr,7*(i-1),0);
Q2ktranfr=[-A1*vz1pr;-BT2'*A1*vz1pr]*f12ktranfr;
QA=Add(QA,Q2ktranfr,7*(j-1),0);    
end   
    
end
k=k+1;
end

end



