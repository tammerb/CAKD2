function [csa,dcsa]=csign(a,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt]=...
    parPart(par);

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


