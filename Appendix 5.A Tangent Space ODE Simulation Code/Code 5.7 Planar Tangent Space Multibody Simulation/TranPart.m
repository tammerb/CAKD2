function [i,j,s1pr,s2pr,v1pr,v2pr]=TranPart(k,PJDT)

i=PJDT(2,k);
j=PJDT(3,k);
s1pr=[PJDT(4,k);PJDT(5,k)];
s2pr=[PJDT(6,k);PJDT(7,k)];
v1pr=[PJDT(9,k);PJDT(10,k)];
v2pr=[PJDT(11,k);PJDT(12,k)];
end
