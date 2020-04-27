function windprofile = windprofileGen(ts,lengthScale_act,overallVariance_act,tf,seed)
% Generates a wind profile under specified conditions
% Inputs :      ts = timestep (in seconds)
%               lengthScale_act = actual length scale (in minutes). SDP
%                                   lookup table should exist for this
%               overallVariance_act = actual overall variance. SDP
%                                       lookup table should exist for this
%               tf = final time (in seconds)
%               seed = positive integer seed for profile randomization
% Outputs:      windprofile = timeseries object where data is the wind
%                               angle, theta, (in degrees). Already has
%                               zero order hold for interpolation

wangChoices = [-45:1/3:45]'; % Predefined wind angles, could change if desired

rng(seed); % for reporducibility
time = [0:ts:tf]'; % time vector
n = size(time,1); %number of samples
ts_min = ts/60; % timescale in minutes

K_act = ts_min/lengthScale_act; 
w_act = sqrt(overallVariance_act*K_act/(.42 + K_act));
markTrans_act = markovSDP(K_act,w_act^2,wangChoices); %Generate Markov matrix

theta_ind = zeros(n,1); % initialize vector
theta_ind(1) = (length(wangChoices)+1)/2; % first angle is always '0' degrees

for i = 1:n-1
    weights = markTrans_act(theta_ind(i),:); % appropriate row of markov matrix
    theta_ind(i+1) = randsample(1:length(wangChoices),1,true,weights); % randomly selects an index based on weights
end
theta = wangChoices(theta_ind); % actual wind angle values

%Create a timeseries object for data
windprofile = timeseries(theta,time,'Name','Theta');
windprofile.DataInfo.Interpolation = 'zoh';
windprofile.DataInfo.Units = 'Degrees';
windprofile.TimeInfo.Units = 'seconds';
windprofile.DataInfo.UserData = ['lt_',num2str(l_t),' ovar_',num2str(ovar),' seed_',num2str(rng)];
end