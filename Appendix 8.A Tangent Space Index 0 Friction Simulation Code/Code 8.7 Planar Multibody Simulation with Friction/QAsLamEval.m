function QAsLam=QAsLamEval(tn,q,qd,Lam,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

QAsLam=zeros(ngc,nc);

k=1;
while k<=nh

if PJDT(1,k)==1     %Revolute joint
[i,j,s1pr,s2pr,R,mus,mud,ms,nm]=RevPart(k,PJDT);
[r1d,ph1d]=qPart(qd,i);
Lamk=[Lam(ms);Lam(ms+1)];
if j==0
[Sfr,Sfrpr]=SfrSfrpr(R*(-ph1d),mus,mud,par);
QAfk1sLam=[zeros(2,2);(R/(norm(Lamk)+10^-6))*Sfr*Lamk'];
QAsLam=Add(QAsLam,QAfk1sLam,3*(i-1),2*(k-1));        
end
if j>=1
[r2d,ph2d]=qPart(qd,j);
[Sfr,Sfrpr]=SfrSfrpr(R*(ph2d-ph1d),mus,mud,par);
QAfk1sLam=[zeros(2,2);(R/(norm(Lamk)+10^-6))*Sfr*Lamk'];
QAsLam=Add(QAsLam,QAfk1sLam,3*(i-1),2*(k-1));
QAfk2sLam=-QAfk1sLam;
QAsLam=Add(QAsLam,QAfk2sLam,3*(j-1),2*(k-1));
end

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
[csLamk1,dcsLamk1]=csign(Lamk(1),par);
f1k=(T1k/d)-d*Lamk(1);
f2k=-T1k/d;
v12=v1pr'*A1'*(-r1d)/d-ph1d*v1pr'*P*s1pr/d;
[Sfr,Sfrpr]=SfrSfrpr(v12,mus,mud,par);
T1ksLamk=-s1pr'*v1pr*[1,0]+[(s2pr'-r1')*A1*v1pr,v1pr'*A1'*v2pr];
f1ksLamk=(1/d)*T1ksLamk-d*[1,0];
f2ksLamk=-(1/d)*T1ksLamk;
[csf1k,dcsf1k]=csign(f1k,par);
[csf2k,dcsf2k]=csign(f2k,par);
absf1ksLamk=(csf1k+f1k*dcsf1k)*f1ksLamk;
absf2ksLamk=(csf2k+f2k*dcsf2k)*f2ksLamk;
QAf1ksLamk=(1/d)*Sfr*[A1*v1pr;v1pr'*P*s1pr]*(absf1ksLamk+absf2ksLamk);
QAsLam=Add(QAsLam,QAf1ksLamk,3*(i-1),2*(k-1));        
end
if j>=1
[r2,ph2]=qPart(q,j);
[r2d,ph2d]=qPart(qd,j);
A2=ATran(ph2);
T1k=-s1pr'*v1pr*Lamk(1)+...
    [(r2'+s2pr'*A2'-r1')*A1*v1pr,v1pr'*A1'*A2*v2pr]*Lamk;
[csLamk1,dcsLamk1]=csign(Lamk(1),par);
f1k=(T1k/d)-d*Lamk(1);
f2k=-T1k/d;
v12=v1pr'*A1'*(r2d+ph2d*P*A2*s2pr-r1d)/d-ph1d*v1pr'*P*s1pr/d;
[Sfr,Sfrpr]=SfrSfrpr(v12,mus,mud,par);
T1ksLamk=-s1pr'*v1pr*[1,0]+[(r2'+s2pr'*A2'-r1')*A1*v1pr,v1pr'*A1'*A2*v2pr];
f1ksLamk=(1/d)*T1ksLamk-d*[1,0];
f2ksLamk=-(1/d)*T1ksLamk;
[csf1k,dcsf1k]=csign(f1k,par);
[csf2k,dcsf2k]=csign(f2k,par);
absf1ksLamk=(csf1k+f1k*dcsf1k)*f1ksLamk;
absf2ksLamk=(csf2k+f2k*dcsf2k)*f2ksLamk;
QAf1ksLamk=(1/d)*Sfr*[A1*v1pr;v1pr'*P*s1pr]*(absf1ksLamk+absf2ksLamk);
QAsLam=Add(QAsLam,QAf1ksLamk,3*(i-1),2*(k-1));
QAf2ksLamk=(1/d)*Sfr*[-A1*v1pr;-v1pr'*A1'*P*A2*s2pr]*...
    (absf1ksLamk+absf2ksLamk);
QAsLam=Add(QAsLam,QAf2ksLamk,3*(j-1),2*(k-1));   
end

end

k=k+1;
end

end

