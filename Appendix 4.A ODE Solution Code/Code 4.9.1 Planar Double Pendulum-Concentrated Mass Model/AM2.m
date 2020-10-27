function M2=AM2(v,mu,par,dat)

[nv,intol,Atol,hmax,hvar]=BparPart(par);
[m1,m2]=AdatPart(dat);
v1=v(1);
v2=v(2);
mu1=mu(1);
mu2=mu(2);

% Enter M2=(M(q,par)mu)sq

M2=sin(v2-v1)*m2*[mu2,-mu2;mu1,-mu1];

end

