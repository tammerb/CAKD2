function [csa,dcsa]=csign(a,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

%Continuous sign function
if abs(a)>vt
    csa=sign(a);
end
    
if abs(a)<=vt
    csa=sin(pi*a/(2*vt));
end

%Derivative of continuous sign function
if abs(a)>vt
    dcsa=0; 
end
    
if abs(a)<=vt
    dcsa=(pi/(2*vt))*cos(pi*a/(2*vt));
end

end





