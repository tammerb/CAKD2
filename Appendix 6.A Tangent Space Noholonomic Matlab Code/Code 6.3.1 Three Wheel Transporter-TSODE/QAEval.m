function QA=QAEval(t,q,qd,par)
%Evaluate Generalized Appoied Force nqx1; example below is for transporter
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par)
phi=q(3);
AT=ATran(phi);
if t<2
    FL=20;
    FR=-20;
end
if t>2
    FL=-20;
    FR=20;
end
if t>4
    FL=0;
    FR=0;
end
Fpr=[0;FL+FR];
QA=[AT*Fpr;FR-FL];



end

