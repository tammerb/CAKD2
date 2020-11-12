function E2=E2eval(t,q,qd,par)
%Enter terms for (Eqd)q*qd, Transporter
[r,phi]=qPart(q);
[rd,phid]=qdPart(qd);
xd=rd(1);
yd=rd(2);
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);
P=[0,-1;1,0];
A=ATran(phi);
uy=[0;1];

E2=-xd*phid*sin(phi)+yd*phid*cos(phi);

end

