function Phiq=PhiqEval(q,par)
%Evaluate Constraint Jacobian

[nq,nh,nv,nu,g,m1,m2,m3,K1,K2,el,mud,mus,vt,...
    utol,Btol,intol,h0]=Partpar(par);

Phiq=[q(1),q(2),0];
    end
