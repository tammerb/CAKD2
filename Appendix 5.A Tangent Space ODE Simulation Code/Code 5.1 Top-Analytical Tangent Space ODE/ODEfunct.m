function vdd=ODEfunct(t,v,vd,par,p0,V,Jpr,m,g,Dcf)

M=AM(v,par,p0,V,Jpr,m,g);
gf=Agf(t,v,vd,par,p0,V,Jpr,m,g,Dcf);
vdd=M\gf;


end

