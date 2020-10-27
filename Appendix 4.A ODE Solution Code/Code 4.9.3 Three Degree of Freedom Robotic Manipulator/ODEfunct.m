function vdd=ODEfunct(t,v,vd,par,dat,integ)

[M,gf,M2,gfsv,gfsvd] = AMg(t,v,vd,vd,par,dat,integ);
vdd=M\gf;


end

