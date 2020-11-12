function csa=csign(a,par)

[nq,nh1,nv1,nu1,nh2,nv2,nu2,g,m1,m2,m3,K1,K2,el,mud,mus,...
    utol,Btol,intol,h]=Partpar(par);

del=h*10;
csa=sign(a);

if (del/2<a)&&(a<del)
    csa=1+0.5*(mus/mud-1)*(1+cos((2*pi/(del))*(a-del/2)));
end
    
if abs(a)<=del/2
    csa=(mus/mud)*sin(pi*a/del);
end

if (-del<a)&&(a<-del/2)
    csa=-1-0.5*((mus/mud)-1)*(1+cos((2*pi/(del))*(a+del/2)));
end



end

