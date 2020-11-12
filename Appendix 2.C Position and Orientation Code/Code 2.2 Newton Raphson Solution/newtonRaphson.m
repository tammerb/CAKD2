function [xFinal] = newtonRaphson(sym_x, xtol, x0, sym_Phi)

xOriginal = x0;
Phiq=jacobian(sym_Phi, sym_x);  %Take Deribative of Phi(x) to obatain Phiq(x)
counter=0;   %Counter for number of iteration
err=xtol+1;  %Set initial error greater than error tolerance

while err > xtol %some thing else Iteration for x, through line 15
    if counter < 20
        Phi_eval=vpa(subs(sym_Phi,sym_x,x0));  %Evaluate F(x) at x0
        Phiq_eval=vpa(subs(Phiq,sym_x,x0));  %Evaluate Fx(x) at x0
        delx = inv(Phiq_eval)*Phi_eval;  %Evaluate delta x
        x0=x0-delx;  %Get a set of new values for iteration/solution
        err=norm(Phi_eval); %Calculate the error using 2 norm
        counter=counter+1;
    else
        disp('The Newton-Raphson has exceeded the max iterations')
        xFinal = xOriginal;
        break
    end
xFinal = x0;    
end
end