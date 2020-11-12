function [xr,xph]=xPart(x,i)

xr=[x(3*(i-1)+1);x(3*(i-1)+2)];
xph=x(3*(i-1)+3);
