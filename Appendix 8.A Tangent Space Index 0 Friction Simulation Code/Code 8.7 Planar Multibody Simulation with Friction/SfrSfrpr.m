function [Sfr,Sfrpr]=SfrSfrpr(v,mus,mud,par)

[nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,...
hvar,NTSDA,NRSDA,vt]=parPart(par);

a=((v/(2*vt))^2+3/4)^2;
b=1-(tanh(4*v/vt))^2;

Sfr=mud*tanh(4*v/vt)+(mus-mud)*(v/vt)/a;   %Friction force function

Sfrpr=(4*mud/vt)*b+(mus-mud)*(1/vt)/a-...
    (mus-mud)*((v^2)/(vt^3))/(a^1.5);       %Derivative of "


end

