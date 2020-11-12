function CtLamsq=CtLamsqeval(t,q,Lam,par)

[r1,phi1,r2,phi2]=qPart(q);
A1=ATran(phi1);
A2=ATran(phi2);
P=[0,-1;1,0];
ipr=[1;0];
jpr=[0;1];
[thet,thetd,thetdd]=Steer(t,par);

Csphi1=[zeros(2,2),3*A1*ipr,zeros(2,2),zeros(2,1);...
    -[-sin(thet),cos(thet)]*A1'*P,0,zeros(1,2),0;...
    -jpr'*A1'*P,0,zeros(1,2),0;zeros(1,6)];
Csphi2=[zeros(2,5),3*A2*ipr;zeros(2,6);zeros(1,2),0,-jpr'*A2'*P,0];

CtLamsq=[zeros(6,2),Csphi1'*Lam,zeros(6,2),Csphi2'*Lam];

end

