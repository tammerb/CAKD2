function [i,j,sipr,sjpr,K,C,el0,F]=STSDATPart(STSDAT,T)

i=STSDAT(1,T);
j=STSDAT(2,T);
sipr=[STSDAT(3,T);STSDAT(4,T);STSDAT(5,T)];
sjpr=[STSDAT(6,T);STSDAT(7,T);STSDAT(8,T)];
K=STSDAT(9,T);
C=STSDAT(10,T);
el0=STSDAT(11,T);
F=STSDAT(12,T);

end

