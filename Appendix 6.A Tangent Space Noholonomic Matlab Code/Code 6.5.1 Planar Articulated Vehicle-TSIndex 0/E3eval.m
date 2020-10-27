function E3=E3eval(t,q,qd,par)

[r1,phi1,r2,phi2]=qPart(q);
[r1d,phi1d,r2d,phi2d]=qdPart(qd);
P=[0,-1;1,0];
A1=ATran(phi1);
A2=ATran(phi2);
jpr=[0;1];
[thet,thetd,thetdd]=Steer(t,par);

E3=[0,0,-phi1d*[-sin(thet),cos(thet)]*A1'*r1d,0,0,0;...
0,0,-phi1d*jpr'*A1'*r1d,0,0,0;0,0,0,0,0,-phi2d*jpr'*A2'*r2d];
end

