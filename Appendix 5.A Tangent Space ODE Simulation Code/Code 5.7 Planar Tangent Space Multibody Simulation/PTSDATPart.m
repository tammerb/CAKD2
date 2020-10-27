function [i,j,sipr,sjpr,K,C,el0,F]=PTSDATPart(PTSDAT,T)

i=PTSDAT(1,T);
j=PTSDAT(2,T);
sipr=[PTSDAT(3,T);PTSDAT(4,T)];
sjpr=[PTSDAT(5,T);PTSDAT(6,T)];
K=PTSDAT(7,T);
C=PTSDAT(8,T);
el0=PTSDAT(9,T);
F=PTSDAT(10,T);

end

