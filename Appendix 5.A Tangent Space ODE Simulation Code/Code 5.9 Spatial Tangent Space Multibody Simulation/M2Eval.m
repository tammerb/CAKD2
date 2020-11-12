function M2=M2Eval(q,mu,SMDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=...
    parPart(par);

M2=zeros(ngc,ngc);

i=1;
while i<=nb
[r,p]=qPart(q,i);
[mur,mup]=qPart(mu,i);
m=SMDT(1,i);
J=diag([SMDT(2,i);SMDT(3,i);SMDT(4,i)]);
G=GEval(p);
Gmu=GEval(mup);
M2i=[zeros(3,7);zeros(4,3),TEval(4*J*G*mup)-4*G'*J*Gmu];
M2=Add(M2,M2i,7*(i-1),7*(i-1));
i=i+1;
end


end

