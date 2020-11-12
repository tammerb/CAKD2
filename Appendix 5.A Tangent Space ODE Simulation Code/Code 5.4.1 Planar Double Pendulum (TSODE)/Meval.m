function M=Meval(q,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

% Enter Mass Matrix M=M(q,par)
I2=eye(2);
M=blkdiag(m1*I2,Jp1,m2*I2,Jp2);

end

