function [r,p,apr,x1,y1]=qPart(q)
%Fill in all q components above and below; example is rolling ellipsoid
r=[q(1);q(2);q(3)];
p=[q(4);q(5);q(6);q(7)];
apr=[q(8);q(9);q(10)];
x1=q(11);
y1=q(12);
%add other components of q 

end

