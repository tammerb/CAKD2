function [B,Biter]= Beval(t,q,B,U,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

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

