function E1=E1eval(t,q,par)
%Enter coefficient matrix of differential constraint; example is transporter 
[r,phi]=qPart(q);
[nq,nh,nd,nv,nu,nw,nx,m,g,FL,FR,integ]=Partpar(par);;    
P=[0,-1;1,0];
E1=[cos(phi),sin(phi),0];



end