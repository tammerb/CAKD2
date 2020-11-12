function Phiq=PhiqEval(q,par,mode)
%Evaluate Constraint Jacobian

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

if mode==1
Phiq=[q(1),q(2),0];
end

if mode==2
Phiq=[q(1),q(2),0;0,1,0];
end

if mode==3
Phiq=[q(1),q(2),0;1,0,0];
end

if mode==4
Phiq=[q(1),q(2),0;0,0,1];
end
    end
