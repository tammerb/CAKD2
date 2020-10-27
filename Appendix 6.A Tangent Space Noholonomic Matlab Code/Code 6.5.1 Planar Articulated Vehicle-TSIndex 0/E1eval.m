function E=E1eval(t,q,par)
 
[r1,phi1,r2,phi2]=qPart(q);
[thet,thetd,thetdd]=Steer(t,par);    
A1=ATran(phi1);
A2=ATran(phi2);
jpr=[0;1];

E=[[-sin(thet),cos(thet)]*A1',3*[-sin(thet),cos(thet)]*jpr,0,0,0;...
jpr'*A1',-3,0,0,0;0,0,0,jpr'*A2',-3];

end

