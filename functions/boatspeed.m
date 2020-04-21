function [u,varargout] = boatspeed(psi,theta,vss_vec)
%BOATSPEED calculates the velocity of the sailboat based on the velocity
%polar
%
%   Assumes a 2 m/s windspeed. Can be
%
%   Inputs:     psi = boat heading angle in degrees.  0 degrees is straight
%                      "up" the course from current position, not 
%                       necessarily pointed at the final waypoint. Positive
%                       is counter-clockwise
%               theta = wind angle in degrees. 0 degrees is wind pointing
%                       from the top of the course downward. CC is positive
%               vss_vec = vector of boatspeeds @ 2 m/s wind from
%                        vss_mat_fine. Keep format of 0:0.1:180 for TWA.
%                         Should be: vss_mat_fine(21,:);
%   Outputs:    u = steady-state longitudinal boat velocity in m/s
%   Optional Output: k = port or starboard tack. k = 1 for a starboard tack
%                        wind over right of boat) and k = 2 for port tack

if nargout>2
    error('Error. Maximum of 2 outputs: boatspeed and tack.')
end 

twa_ref = 0:0.1:180;
twa_act = psi - theta;
u = interp1(twa_ref,vss_vec,abs(twa_act));
if nargout==2 % For tack of path
    varargout{1} = 1 + (psi>theta); %
end

end
