function [B,Biter]= Beval(t,q,B,U,par)

[nq,nh,utol,Btol,intol,Atol,m1,m2,Jp1,Jp2,g,...
    K1,K2,K3,C1,C2,n1,n2]=parPart(par);

i=1;
Cnorm=Btol+1;
I=eye(nh);
Phisq=P1(t,q,par);
while Cnorm>Btol
i=i+1;
delB=-B*Phisq*U*B+B;
B=B+delB;
C=B*Phisq*U-I;
Cnorm=norm(C);
end
Biter=i-1;

end

