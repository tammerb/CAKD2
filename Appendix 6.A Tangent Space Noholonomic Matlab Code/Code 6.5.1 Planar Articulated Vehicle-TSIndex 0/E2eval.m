function E2=E2eval(t,q,x,par)

[r1,phi1,r2,phi2]=qPart(q);
[xr1,xphi1,xr2,xphi2]=qPart(x);
P=[0,-1;1,0];
A1=ATran(phi1);
A2=ATran(phi2);
jpr=[0;1];
[thet,thetd,thetdd]=Steer(t,par);

E2=[0,0,-[-sin(thet),cos(thet)]*A1'*P*xr1,0,0,0;0,0,-jpr'*A1'*P*xr1,0,0,0;...
    0,0,0,0,0,-jpr'*A2'*P*xr2];

end

