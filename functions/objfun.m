function J = objfun(u,x,y,theta,xd,k,sdp_ctg,vss_vec,ts,gate_w)
%OBJFUN Calculates the objective function for a given control sequence
%   Objective function derived for minimizing time between the current
%   position and a specified waypoint
%   Inputs: u   = control sequence
%           x   = current x-position
%           y   = current y-position
%           theta = current wind angle
%           xd  = x-position of desired waypoint (from SDP lookup table)
%           k   = previous tack. Where 1 = starboard tack and 2 = port tack
%           sdp_ctg = 2-element vector from ctg(xd_ind,theta_ind,stage#,:)
%                     representing the cost to go from the destination when
%                     arriving on either a starboard or port tack. xd_ind
%                     and theta_ind are the indices of the values xd and
%                     theta from their respective choices. For example, xd
%                     options from sdp are -800:2:800, so xd_ind(xd = 0) =
%                     401, or theta_ind(theta = -35) = 3; If theta is not
%                     contained within -45:5:45, use the index of the
%                     corresponding bin where bin EDGES are at -47.5:5:47.5
%                     First element should correspond to the starboard tack
%                     and the second element from the port tack
%           vss_vec = vector of boatspeeds @ 2 m/s wind from
%                        vss_mat_fine. Keep format of 0:0.1:180 for TWA.
%                         Should be: vss_mat_fine(21,:);
%           ts = timestep
%           gate_w = total gatewidth (half on either side) of next waypoint 
%                    at boundary between SDP stages (in meters)
%
%   Output: J = objective function value


% Boat velocity and position
yDisc = 10; % stage length
tPen = 5; % tacking penalty

yd = ceil(y/yDisc)*yDisc; % next y-waypoint
[vel,tack] = boatspeed(u,theta,vss_vec); % total velocity and tack at every timestep
yVel = vel.*cosd(u);    % y velocity
yProg = yVel*ts;        % y progress made during each timestep
yPos = y + cumsum(yProg);       % y position after each timestep
ts_end = find(yPos>=yd,1);      % crosses into next stage after this many timesteps
if isempty(ts_end)
    ts_end = length(u);
    no_finish = true; % used later for calculating penalty from finish
else
    no_finish = false;
end
xVel = vel(1:ts_end).*sind(u(1:ts_end)); % x velocity
xProg = xVel*ts; % x progress
xPos = x + cumsum(xProg); % x position until timestep where boat crosses into next stage

% Where position where boat crosses into next stage
ts_end_per = (yd - yPos(ts_end-1))/yProg(ts_end); % proportion of time on last stage to enter into next stage
xf = xPos(ts_end - 1) + xProg(ts_end)*ts_end_per; % x position as boat crosses into next stage

% Time to reach next SDP stage
tf = (ts_end - 1 + ts_end_per)*ts; % elapsed time to next stage

% Calculate distance for penalty at finish
if (~no_finish && abs(xf-xd)<=gate_w/2) % finished within gate
    d = 0;
elseif (~no_finish && abs(xf-xd)>=gate_w/2) % finished outside of gate
    d = min(abs(xf-(xd - gate_w/2)),abs(xf - (xd + gate_w/2)));
else
    d = norm([xf yPos(end)] - [xd yd]);
end

% Calculate total number of tacking penalties
tack_aug = [k; tack(:)];
numtacks = nnz(diff(tack_aug));

% Objective function value
J = tf + numtacks*tPen + sdp_ctg(tack(end)) + distpen(d);

end

function y = distpen(d)
% Calculates the distance penalty
a = 10; % coeff on quadratic distance penalty
y = a*d^2;
end




