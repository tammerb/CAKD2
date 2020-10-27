function D=Deval(t,q,U,V,B,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

Phisq=P1(t,q,par);
D=(eye(nq)-U*B*Phisq)*V;
end

