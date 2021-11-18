function [Wf] = SRP_eclipse_model(alpha, D)
% This funtction gives the value of the Weigth factor to use for the solar
% radiation when the spacecraft is located in the umbra or penumbra zone
% during eclipses.
%
% The input are:
%
% [alpha]: it's the angle beteween the Sun-Earth direction and the
% spacecraft position which is found using the vector which connects the SC
% and the center  of the Earth; must be in [rad], I don't know why...
% 
% [D]: it' the distance beteween the Earth center and the projection of the
% position of the SC on the direction Sun-Earth
%
%
% This model works for the Sun-Earth couple eclipse, but can works also for
% other celestial body changing the input paramteres and giving the
% characteristic of the planet

    Wf = []; %[-]
    Re = 6371; % earth radius[km]
    Rs = 696340; % sun radius [km]
    Ds_e = 150000000; % mean sun-earth distance [km]
    Ys = D*tan(alpha); % S/C-(sun_earth direction) distance [km]
    Ys_abs = abs(Ys);
    De = 2*Re; Dp = 2*Re; %[km]
    Ds = 2*Rs;
    
%% umbral geometry
    X_u = (Dp*Ds_e)/(Ds-Dp); %[km]
    alpha_u =abs( asin((Dp/(2*X_u))));%[rad]
%% penumbral geometry
    X_p = (Dp*Ds_e)/(Ds+Dp);%[km]
    alpha_p = abs(asin((Dp/(2*X_p))));%[rad]
%% shadow terminator
    Umbra_terminator = ((X_u-D)*tan(alpha_u));%[km]
    Penumbra_terminator = ((X_p+D)*tan(alpha_p));%[km]
    
if Ys_abs >= (Penumbra_terminator)
    
    Wf = 1;
 
else if Ys_abs > (Umbra_terminator) && Ys_abs < (Penumbra_terminator)
    
        Wf = ( (Ys_abs-(Umbra_terminator))/((Penumbra_terminator)-(Umbra_terminator)) );

else if Ys_abs <= (Umbra_terminator)
        
        Wf = 0;
 
        end
        end
end


end