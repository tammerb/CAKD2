function [u,Iteru] = usolv(u,v,q0,V,U,B,par,mode,q1b,q2b,q3b)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

err=utol+1;
i=1;
q=q0+V*v-U*u;
while err>utol
delu=B*PhiEval(q,par,mode,q1b,q2b,q3b);
u=u+delu;
q=q0+V*v-U*u;
i=i+1;
err=norm(PhiEval(q,par,mode,q1b,q2b,q3b));
end
Iteru=i-1;
end

