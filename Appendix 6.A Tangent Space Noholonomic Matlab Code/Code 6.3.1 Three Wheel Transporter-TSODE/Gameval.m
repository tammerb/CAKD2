function Gam=Gameval(t,q,qd,par)
%Enter expressions for eqqd, et, Ptq, Ptt, and Gam; planar art veh
[r,phi]=qPart(q)
[rd,phid]=qdPart(qd)
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);

E2=E2eval(t,q,qd,par);

Gam=-E2;



end



