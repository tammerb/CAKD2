function QA=QAeval(t,q,qd,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

uy=[0;1;0];

QA=[-m*g*uy;zeros(4,1)];
end

