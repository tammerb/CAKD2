function [i,j,s1pr,s2pr,d,ms,nm]=DistPart(k,SJDT)

i=SJDT(2,k);
j=SJDT(3,k);
s1pr=[SJDT(4,k);SJDT(5,k);SJDT(6,k)];
s2pr=[SJDT(7,k);SJDT(8,k);SJDT(9,k)];
d=SJDT(10,k);
ms=SJDT(25,k);
nm=SJDT(26,k);
end

