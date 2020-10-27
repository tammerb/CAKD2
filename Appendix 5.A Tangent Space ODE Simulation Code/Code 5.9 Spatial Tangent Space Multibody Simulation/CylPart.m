function [i,j,s1pr,s2pr,ux1pr,uz1pr,ux2pr,uz2pr]=CylPart(k,SJDT)

i=SJDT(2,k);
j=SJDT(3,k);
s1pr=[SJDT(4,k);SJDT(5,k);SJDT(6,k)];
s2pr=[SJDT(7,k);SJDT(8,k);SJDT(9,k)];
ux1pr=[SJDT(11,k);SJDT(12,k);SJDT(13,k)];
uz1pr=[SJDT(14,k);SJDT(15,k);SJDT(16,k)];
ux2pr=[SJDT(17,k);SJDT(18,k);SJDT(19,k)];
uz2pr=[SJDT(20,k);SJDT(21,k);SJDT(22,k)];

end

