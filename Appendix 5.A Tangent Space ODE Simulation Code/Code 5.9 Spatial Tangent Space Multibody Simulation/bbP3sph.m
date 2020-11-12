function [P31,P32]=bbP3sph(tn,q,qd,par)


[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

%P3 spherical contribution is zero, no nonzero terms to add
P31=zeros(1,7);
P32=zeros(1,7);
end









