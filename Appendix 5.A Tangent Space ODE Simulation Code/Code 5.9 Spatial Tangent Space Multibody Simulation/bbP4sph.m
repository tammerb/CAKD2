function [P411,P412,P422]=bbP4sph(i,j,s1pr,s2pr,tn,q,etak,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,h0,hvar,NTSDA]=parPart(par);

I3=eye(3);

%P412 is identically zero
P412=zeros(7,7);

if j==0
P411=[zeros(3,7);zeros(4,3),-KEval(s1pr,etak)];
P422=zeros(7,7);
end

if j>=1
P411=[zeros(3,7);zeros(4,3),-KEval(s1pr,etak)];
P422=[zeros(3,7);zeros(4,3),KEval(s2pr,etak)];
end

end









