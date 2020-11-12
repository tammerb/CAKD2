function [Ssq,Ssqd]=Ssqqdeval(q,qd,par)

[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);
[r,p,a]=qPart(q);
[rd,pd,ad]=qdPart(qd);

Gbarpd=Gbareval(pd);
Gbarp=Gbareval(p);

Ssq=[zeros(3,9);zeros(4,3),8*Gbarpd'*Jpr*Gbarpd,zeros(4,2);zeros(2,9)];
Ssqd=[zeros(3,9);zeros(4,3),8*(Tbareval(Jpr*Gbarpd*p)-Gbarpd'*Jpr*Gbarp),...
    zeros(4,2);zeros(2,9)];

end

