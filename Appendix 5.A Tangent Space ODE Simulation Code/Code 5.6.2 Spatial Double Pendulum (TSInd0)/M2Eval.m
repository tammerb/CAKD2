function M2=M2Eval(q,b,par)

p1=[q(1);q(2);q(3);q(4)];
p2=[q(5);q(6);q(7);q(8)];
Gp1=GEval(p1);
Gp2=GEval(p2);
bp1=[b(1);b(2);b(3);b(4)];
bp2=[b(5);b(6);b(7);b(8)];
Gbp1=GEval(bp1);
Gbp2=GEval(bp2);
T1=TEval(Gp1*bp1);
T2=TEval(Gp2*bp2);
M2=(8*par(3)/5)*blkdiag(-Gp1'*Gbp1+T1,-Gp2'*Gbp2+T2,zeros(3,3));


end

