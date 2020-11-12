function [B,Biter]=CorrectB(tn,q,B,U,SJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

i=1;
Cnorm=Btol+1;
I=eye(nu);
Phiq=PhiqEval(tn,q,SJDT,par);
while Cnorm>Btol
delB=-B*Phiq*U*B+B;
B=B+delB;
C=B*Phiq*U-I;
Cnorm=norm(C);
i=i+1;
end
Biter=i-1;
end

