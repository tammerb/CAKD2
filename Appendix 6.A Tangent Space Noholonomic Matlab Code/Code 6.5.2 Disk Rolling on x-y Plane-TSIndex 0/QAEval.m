function QA=QAEval(t,q,qd,par)
%Evaluate Generalized Appoied Force nqx1; example below is for rolling ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,Jpr,utol,Btol,Htol,intol,h]=Partpar(par);

uz=[0;0;1];

QA=[-m*g*uz;zeros(6,1)];
end

