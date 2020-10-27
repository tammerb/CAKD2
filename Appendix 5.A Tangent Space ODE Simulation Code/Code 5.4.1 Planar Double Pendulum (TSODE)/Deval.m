function D=Deval(t,q,U,V,B,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

Phisq=P1(t,q,par);
D=(eye(nq)-U*B*Phisq)*V;
end

