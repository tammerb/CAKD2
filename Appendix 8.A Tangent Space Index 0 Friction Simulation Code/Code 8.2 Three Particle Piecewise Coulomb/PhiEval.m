function Phi=PhiEval(q,par,mode,q1b,q2b,q3b)
%Evaluate constraint function

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);
if mode==1
Phi=(q(1)^2+q(2)^2-el^2)/2;
end

if mode==2
Phi=[(q(1)^2+q(2)^2-el^2)/2;q(2)-q2b];
end

if mode==3
Phi=[(q(1)^2+q(2)^2-el^2)/2;q(1)-q1b];
end

if mode==4
Phi=[(q(1)^2+q(2)^2-el^2)/2;q(3)-q3b];
end

end

