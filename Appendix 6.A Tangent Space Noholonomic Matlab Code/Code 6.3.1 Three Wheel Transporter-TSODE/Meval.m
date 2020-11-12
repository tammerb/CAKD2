function M=Meval(q,par)
%Enter Mass Matrix nqxnq; example below is transporter
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);
P=[0,-1;1,0];
phi=q(3);
AT=ATran(phi);
spc=[0;1];
M=[m*eye(2),m*P*AT*spc;(m*P*AT*spc)',m];


end

