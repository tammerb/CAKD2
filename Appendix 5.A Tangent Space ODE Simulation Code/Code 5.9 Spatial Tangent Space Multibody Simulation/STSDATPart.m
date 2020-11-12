function [i,j,sipr,sjpr,K,C,el0,F]=STSDATPart(STSDAT,T)

a=STSDAT(:,T);
i=a(1);
j=a(2);
sipr=[a(3);a(4);a(5)];
sjpr=[a(6);a(7);a(8)];
K=a(9);
C=a(10);
el0=a(11);
F=a(12);

end

