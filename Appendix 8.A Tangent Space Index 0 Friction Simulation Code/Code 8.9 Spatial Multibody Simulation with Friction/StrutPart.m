function [i,j,s1pr,s2pr,ux2pr,uy2pr]=StrutPart(k,SJDT)

i=SJDT(2,k);
j=SJDT(3,k);
s1pr=[SJDT(4,k);SJDT(5,k);SJDT(6,k)];
s2pr=[SJDT(7,k);SJDT(8,k);SJDT(9,k)];
ux2pr=[SJDT(17,k);SJDT(18,k);SJDT(19,k)];
uz2pr=[SJDT(20,k);SJDT(21,k);SJDT(22,k)];
uy2pr=atil(uz2pr)*ux2pr;
end

