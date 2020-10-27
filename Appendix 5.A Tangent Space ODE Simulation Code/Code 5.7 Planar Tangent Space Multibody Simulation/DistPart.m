function [i,j,s1pr,s2pr,d]=DistPart(k,PJDT)

i=PJDT(2,k);
j=PJDT(3,k);
s1pr=[PJDT(4,k);PJDT(5,k)];
s2pr=[PJDT(6,k);PJDT(7,k)];
d=PJDT(8,k);
end

