function [i,j,s1pr,s2pr,vz1pr,vz2pr]=UnivPart(k,SJDT)

i=SJDT(2,k);
j=SJDT(3,k);
s1pr=[SJDT(4,k);SJDT(5,k);SJDT(6,k)];
s2pr=[SJDT(7,k);SJDT(8,k);SJDT(9,k)];
vz1pr=[SJDT(14,k);SJDT(15,k);SJDT(16,k)];
vz2pr=[SJDT(20,k);SJDT(21,k);SJDT(22,k)];
end

