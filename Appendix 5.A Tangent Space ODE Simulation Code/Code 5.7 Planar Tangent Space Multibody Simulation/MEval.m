function M=MEval(PMDT,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

M=zeros(ngc,ngc);
I2=eye(2);
i=1;
while i<=nb
m=PMDT(1,i);
J=PMDT(2,i);
Mi=blkdiag(m*I2,J);
M=Add(M,Mi,3*(i-1),3*(i-1));
i=i+1;
end

end

