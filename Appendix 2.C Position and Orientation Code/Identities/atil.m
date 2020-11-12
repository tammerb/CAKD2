function atil=atil(a)
% The tilda operator from eq 2.1.22.
% a = [a_x,a_y,a_z]'
% atil(a) = [0 -a_z a_y]
%           [a_z 0 -a_x]
%           [-a_y a_x 0]
atil=[0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
end
