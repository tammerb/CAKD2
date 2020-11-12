function dcsa=dcsign(a,par)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h,mode,q1b,q2b]=Partpar(par);

dcsa=0;

if (delta/2<a)&&(a<delta)
    dcsa=-0.5*(mus/mud-1)*((2*pi/(delta))*sin((2*pi/(delta))*(a-delta/2))); 
end
    
if abs(a)<=delta/2
    dcsa=(mus*pi/(mud*delta))*cos(pi*a/delta);
end

if (-delta<a)&&(a<-delta/2)
    dcsa=0.5*((mus/mud)-1)*((2*pi/(delta))*sin((2*pi/(delta))*(a+delta/2)));
end



end

