function PqTLsq=PhiqTLamsq(Lam,par,mode)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

if mode==1
PqTLsq=[Lam,0,0;0,Lam,0;0,0,0];
end
if mode>1
PqTLsq=[Lam(1),0,0;0,Lam(1),0;0,0,0];
end

end

