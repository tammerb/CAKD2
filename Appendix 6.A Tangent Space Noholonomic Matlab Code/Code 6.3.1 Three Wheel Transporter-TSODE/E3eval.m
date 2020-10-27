function E3=E3eval(t,q,qd,par)
[r,phi]=qPart(q)
[rd,phid]=qdPart(qd)
[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ]=Partpar(par);
P=[0,-1;1,0];
A=ATran(phi);

uy=[0;1];
ux=[1;0];
[thet,thetd,thetdd]=Steer(t,par)

E3=[0;thetd*[-sin(thet),cos(thet)]*(P'*A'*rd+phid*ux)];
end

