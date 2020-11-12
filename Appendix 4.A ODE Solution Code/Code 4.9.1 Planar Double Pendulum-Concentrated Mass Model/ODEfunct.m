function vdd=ODEfunct(t,v,vd,par,dat)

M=AM(v,par,dat);
gf=Agf(t,v,vd,par,dat);
vdd=M\gf;


end

