function QA=QAEval(t,q,qd,par)
%Evaluate Generalized Appoied Force nqx1; example below is for rolling ellipsoid
[nq,nh,nd,nv,nu,nw,nx,m,g,F,K,eps,amp,om]=parPart(par);
r=[q(1);q(2);q(3)];
uz=[0;0;1];
f=[F*sin(t);F*cos(2*t);-m*g];
QA=[f-K*(r);0;0;0;0;0;0;0;0;0];


end

