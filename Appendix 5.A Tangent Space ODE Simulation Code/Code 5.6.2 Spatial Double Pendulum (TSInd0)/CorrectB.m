function [B,Biter]=CorrectB(q,B,U,Btol)
i=1;
Cnorm=Btol+1;
I=eye(3);
Phiq=PhiqEval(q);
while Cnorm>Btol
i=i+1;
delB=-B*Phiq*U*B+B;
B=B+delB;
C=B*Phiq*U-I;
Cnorm=norm(C);
Biter=i-1;
end
end

