function [B,Biter]=BEval(tn,q,B,U,PJDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

i=1;
Cnorm=Btol+1;
I=eye(nu);
Phiq=PhiqEval(tn,q,PJDT,par);
while Cnorm>Btol
delB=-B*Phiq*U*B+B;
B=B+delB;
C=B*Phiq*U-I;
Cnorm=norm(C);
i=i+1;
end
Biter=i-1;
end

