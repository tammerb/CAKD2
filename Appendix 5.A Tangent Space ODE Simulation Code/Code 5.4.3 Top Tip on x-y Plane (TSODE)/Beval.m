function [B,Biter]= Beval(t,q,B,U,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,K,C]=parPart(par);

i=1;
%Changed C to C1 to avoid conflict with damping coefficient
C1norm=Btol+1;
I=eye(nh);
Phisq=P1(t,q,par);
while C1norm>Btol
i=i+1;
delB=-B*Phisq*U*B+B;
B=B+delB;
C1=B*Phisq*U-I;
C1norm=norm(C1);
end
Biter=i-1;

end

