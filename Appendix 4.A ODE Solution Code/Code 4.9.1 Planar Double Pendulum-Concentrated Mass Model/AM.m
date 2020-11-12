function M=AM(v,par,dat)

[nv,intol,Atol,hmax,hvar]=BparPart(par);
[m1,m2]=AdatPart(dat);
v1=v(1);
v2=v(2);

% Enter Mass Matrix M=M(v,par)

M=[m1+m2,m2*cos(v2-v1);m2*cos(v2-v1),m2];

end

