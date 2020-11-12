function Gam=Gameval(t,q,qd,par)


P2=P2eval(t,q,qd,par);
E2=E2eval(t,q,qd,par);

[r1,phi1,r2,phi2]=qPart(q);
[thet,thetd,thetdd]=Steer(t,par);    
A1=ATran(phi1);
A2=ATran(phi2);
jpr=[0;1];
c=cos(thet);
s=sin(thet);
Est=-thetd*[[c,s]*A1',3*[c,s]*jpr,0,0,0;zeros(2,6)];

Gam=[P2*qd;E2*qd+Est*qd];

end



