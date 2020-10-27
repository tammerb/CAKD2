function [Gamsq,Gamsqd] = Gamsqqdeval(t,q,qd,par)

P2=P2eval(t,q,qd,par);
E2=E2eval(t,q,qd,par);
E3=E3eval(t,q,qd,par);
P3=P3eval(t,q,qd,par);

[r1,phi1,r2,phi2]=qPart(q);
[r1d,phi1d,r2d,phi2d]=qdPart(qd);
[thet,thetd,thetdd]=Steer(t,par);    
A1=ATran(phi1);
A2=ATran(phi2);
jpr=[0;1];
c=cos(thet);
s=sin(thet);
Est=-thetd*[[c,s]*A1',3*[c,s]*jpr,0,0,0;zeros(2,6)];

P=[0,-1;1,0];
Estqdsq=thetd*[0,0,[c,s]*A1'*P*r1d,0,0,0;zeros(2,6)];

Gamsq=[P3;E3+Estqdsq];
Gamsqd=[2*P2;2*E2+Est];


end

