function [i,j,vx1pr,vy1pr,vx2pr]=RotDrPart(k,SJDT)

%Data derived from corresponding revolute or cylindriocal joint data

i=SJDT(2,k);
j=SJDT(3,k);
vx1pr=[SJDT(11,k);SJDT(12,k);SJDT(13,k)];
vz1pr=[SJDT(14,k);SJDT(15,k);SJDT(16,k)];
vy1pr=atil(vz1pr)*vx1pr;
vx2pr=[SJDT(17,k);SJDT(18,k);SJDT(19,k)];


end
