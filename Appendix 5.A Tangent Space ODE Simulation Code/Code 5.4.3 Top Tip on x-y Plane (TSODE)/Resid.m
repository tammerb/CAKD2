function R=Resid(tn,ue,v,vd,vdd,q0,V,U,B,par)

[u,Iteru]=ueval(tn,ue,v,q0,V,U,B,par);
q=q0+V*v-U*u;
[B,Biter]=Beval(tn,q,B,U,par);
D=Deval(tn,q,U,V,B,par);
[Pst,Pstt,Pstq,Psttq]=P5(tn,q,par);
qd=D*vd-U*B*Pst;
QA=QAeval(tn,q,qd,par);
S=Seval(q,qd,par);
M=Meval(q,par);
Gamma=Gameval(tn,q,qd,par);
R=D'*M*D*vdd-D'*(M*U*B*Gamma+S+QA);


end

