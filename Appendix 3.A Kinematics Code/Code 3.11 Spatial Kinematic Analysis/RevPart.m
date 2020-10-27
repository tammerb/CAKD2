function [i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr]=RevPart(k,SJDT)

i=SJDT(2,k);
j=SJDT(3,k);
s1pr=[SJDT(4,k);SJDT(5,k);SJDT(6,k)];
s2pr=[SJDT(7,k);SJDT(8,k);SJDT(9,k)];
vx1pr=[SJDT(11,k);SJDT(12,k);SJDT(13,k)];
vz1pr=[SJDT(14,k);SJDT(15,k);SJDT(16,k)];
vx2pr=[SJDT(17,k);SJDT(18,k);SJDT(19,k)];
vz2pr=[SJDT(20,k);SJDT(21,k);SJDT(22,k)];

end



