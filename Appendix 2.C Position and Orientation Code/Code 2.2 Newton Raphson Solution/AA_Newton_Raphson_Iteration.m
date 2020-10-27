% Newton Raphson Solution of k Equaitons Phi(q)=0 for k Variables q

qtol=0.001;  % Error tolerance in satisfying Phi(q)=0
q=q0;   % Initial solution estimate

i=1;    %Iteration counter
err=qtol+1;
while err >qtol  %Iteration for q, through line 15
Phi=PhiEval(q);
Phiq=PhiqEval(q);
delq=-Phiq\(Phi);  %Newton-Raphson iteration
q=q+delq;
err=norm(Phi);
i=i+1;    
end
iter=i-1;  %Report number of iteations