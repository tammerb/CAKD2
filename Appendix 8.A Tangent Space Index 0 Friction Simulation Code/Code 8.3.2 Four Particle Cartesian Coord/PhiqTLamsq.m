function PqTLsq = PhiqTLamsq(Lam)

PqTLsq=[Lam(7),0,0,-Lam(7),zeros(1,6);0,Lam(7),0,0,-Lam(7),zeros(1,5);...
    0,0,Lam(7),0,0,-Lam(7),zeros(1,4);...
    -Lam(7),0,0,Lam(7)+Lam(8),0,0,-Lam(8),0,0,0;...
    0,-Lam(7),0,0,Lam(7)+Lam(8),0,0,-Lam(8),0,0;
    0,0,-Lam(7),0,0,Lam(7)+Lam(8),0,0,-Lam(8),0;...
    0,0,0,-Lam(8),0,0,Lam(8),0,0,0;
    0,0,0,0,-Lam(8),0,0,Lam(8),0,0;...
    0,0,0,0,0,-Lam(8),0,0,Lam(8),0;zeros(1,10)];


end
