function QA=QAeval(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

uz=[0;0;1];

QA=[-m*g*uz;zeros(4,1)];
end

