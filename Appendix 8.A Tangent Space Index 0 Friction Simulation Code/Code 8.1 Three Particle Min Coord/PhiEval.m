function Phi = PhiEval(q,par)
%Evaluate constraint function

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

Phi=(q(1)^2+q(2)^2-el^2)/2;
end

