% Script for looping simulations based on wind parameters and/or rng seed
close; clear;

%% PAY ATTENTION TO COMMENT YOU MADE AT THE BOTTOM OF THE SCRIPT!!!
tic
% Specify wind parameters
lt = 20; % length scale (minutes)
ovar = 1000; % overall variance (deg^2)
% Specify rng seeds
rng_seeds = 1:4;

MPCtimes = zeros(length(rng_seeds),1); % init
MPCtacktotals = MPCtimes;
SDPtimes = zeros(length(rng_seeds),1); % init
SDPtacktotals = SDPtimes;
tscStore = cell(2,length(rng_seeds)); % store tsccollections
tscMPCStore = cell(1,length(rng_seeds)); % store tsccollections
tscSDPStore = cell(1,length(rng_seeds)); % store tsccollections
parfor i = 1:length(rng_seeds)
%      fprintf('rng = %d \n',i)
    tscMPC = MPCsim(lt,ovar,rng_seeds(i));
    tscSDP = SDPsim(lt,ovar,rng_seeds(i));
    MPCtimes(i) = tscMPC.Time(end);
    MPCtacktotals(i) = str2double(tscMPC.TimeInfo.UserData);
    SDPtimes(i) = tscSDP.Time(end);
    SDPtacktotals(i) = str2double(tscSDP.TimeInfo.UserData);
    tscStore(:,i) = {tscMPC;tscSDP};
end
toc

% Add the appropriate lines here to automatically save results and append a
% results mat file just in case something crashes. 
% Ex:     save('testsave.mat','tscStore_recentruns','-append')