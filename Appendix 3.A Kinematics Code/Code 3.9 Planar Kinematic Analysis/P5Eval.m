function [Pf,Pfst,Pfstt]=P5Eval(tn,par)

[nb,ngc,nh,nhc,nd,qtol,app]=parPart(par);

% Enter Constraint t derivatives of P(t,par);

if app==1   %Slider-Crank
omega=1;
PDf=[omega*tn];        %Enter nd driver functions of tn
Pf=[zeros(nhc,1);PDf];  %Time dependent driver vector
PDfd=[omega];       %Enter nd first derivatives of driver functions of tn
Pfst=[zeros(nhc,1);PDfd];  %Time dependent driver velocity vector
PDfdd=[0];      %Enter nd second derivatives of driver functions of tn   
Pfstt=[zeros(nhc,1);PDfdd];  %Time dependent driver acceleration vector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if app==2   %Quick Return
omega=4;
PDf=[omega*tn];        %Enter driver functions of tn
Pf=[zeros(nhc,1);PDf];
PDfd=[omega];       %Enter first derivatives of driver functions of tn
Pfst=[zeros(nhc,1);PDfd];
PDfdd=[0];      %Enter second derivatives of driver functions of tn   
Pfstt=[zeros(nhc,1);PDfdd];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if app==3   %Windshield Wiper
omega=2;
PDf=[omega*tn];        %Enter driver functions of tn
Pf=[zeros(nhc,1);PDf];
PDfd=[omega];       %Enter first derivatives of driver functions of tn
Pfst=[zeros(nhc,1);PDfd];
PDfdd=[0];      %Enter second derivatives of driver functions of tn   
Pfstt=[zeros(nhc,1);PDfdd];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
