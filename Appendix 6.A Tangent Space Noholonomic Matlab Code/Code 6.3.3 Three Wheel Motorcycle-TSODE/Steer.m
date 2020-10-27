function [thet,thetd,thetdd]=Steer(t,par)
%Defines steer angle for selected mode
[nq,nh,nd,nv,nu,nw,nx,m,g,F,k,c,phi,amp,om,sf,utol,Btol,Htol,...
    mode,integ,hmax]=parPart(par);
thet=amp*sin(om*t);
thetd=amp*om*cos(om*t);
thetdd=-amp*(om^2)*sin(om*t);
    
    if mode==1
        if om*t>2*pi
            thet=0
            thetd=0
            thetdd=0;
        end
    end
        
     if mode==2
         if om*t>pi/2
            thet=amp;
            thetd=0;
            thetdd=0;
         end
     end
            


end

