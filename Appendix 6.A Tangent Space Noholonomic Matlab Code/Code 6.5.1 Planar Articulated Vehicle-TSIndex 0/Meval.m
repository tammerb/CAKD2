function M=Meval(q,par)
%Enter Mass Matrix nqxnq; example below is for planar articulated vehicle
[nq,nh,nd,nv,nu,nw,nx,m,g,amp,om,mode,integ]=Partpar(par);

M=m*eye(6);


end

