function [i,j,sipr,sjpr,K,C,el0,F]=STSDATPart(STSDAT,T)

i=STSDAT(1);
j=STSDAT(2);
sipr=[STSDAT(3);STSDAT(4);STSDAT(5)];
sjpr=[STSDAT(6);STSDAT(7);STSDAT(8)];
K=STSDAT(9);
C=STSDAT(10);
el0=STSDAT(11);
F=STSDAT(12);

end

