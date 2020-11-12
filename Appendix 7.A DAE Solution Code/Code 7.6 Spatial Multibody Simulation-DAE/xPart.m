function [xr,xp]=xPart(x,i)

xr=[x(7*(i-1)+1);x(7*(i-1)+2);x(7*(i-1)+3)];
xp=[x(7*(i-1)+4);x(7*(i-1)+5);x(7*(i-1)+6);x(7*(i-1)+7)];

end
