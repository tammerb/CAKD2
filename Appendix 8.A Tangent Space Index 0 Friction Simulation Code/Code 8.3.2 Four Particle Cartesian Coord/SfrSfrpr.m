function [Sfr,Sfrpr]=SfrSfrpr(v,par)

[nq,nh,nv,nu,g,m1,m2,m3,m4,K1,K2,K3,K4,el1,el2,mud,mus,vt,...
    utol,Btol,intol,h0,FR]=Partpar(par);

a=((v/(2*vt))^2+3/4)^2;
b=1-(tanh(4*v/vt))^2;

Sfr=mud*tanh(4*v/vt)+(mus-mud)*(v/vt)/a;   %Friction force function

Sfrpr=(4*mud/vt)*b+(mus-mud)*(1/vt)/a-...
    (mus-mud)*((v^2)/(vt^3))/(a^1.5);       %Derivative of "


end

