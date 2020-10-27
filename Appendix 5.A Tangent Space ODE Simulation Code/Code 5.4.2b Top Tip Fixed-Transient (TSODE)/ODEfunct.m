function [f,E,B]=ODEfunct(t,v,vd,ue,U,V,B,q0,par)

[nq,nh,utol,Btol,intol,Atol,m,g,Ixy,Iz,Dcf]=parPart(par);

[u,Iteru]=ueval(t,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(t,q,B,U,par);
D=Deval(t,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(t,q,par);
qd=D*vd-U*B*Pst;
M=Meval(q,par);
E=D'*M*D;
QA=QAeval(t,q,qd,par);
S=Seval(q,qd,par);
RHS=D'*(M*U*B*Gameval(t,q,qd,par)+S+QA);
f=E\RHS;


end

