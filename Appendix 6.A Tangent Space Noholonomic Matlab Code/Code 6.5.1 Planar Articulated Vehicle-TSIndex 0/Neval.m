function N=Neval(t,q,x,par)

P2=P2eval(t,q,x,par);
E2=E2eval(t,q,x,par);

N=[P2;E2];


end

