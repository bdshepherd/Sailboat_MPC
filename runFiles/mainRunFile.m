clear
clc
% main run file

%% generate stochastic wind
% length scale of stochastic wind
lengthScale_act = 10;
% variance of stochastic wind
overallVariance_act = 500;
% time step for wind as well as MPC
timeStep = 0.5;
% final time for wind gereration as wind it is generated as a timeseries
tFinal = 60*60; % minutes*60
% random number generator seed for stochastic wind generation
rngSeed = 1;
% run the stochastic wind generation function
windprofile = windprofileGen(timeStep,lengthScale_act,overallVariance_act,...
    tFinal,rngSeed);

%% load SDP 'cost to go' and 'optimal path' matrices and vss_mat
% load the appropriate mat file
load(['sdp_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act),'.mat']);
% evaluate the cost to go matrix
ctg = eval(['ctg_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act)]);
% evaluate the optimal path matrix
optimalPath = eval(['optpath_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act)]);
% load the vss_mat matrix
load('vss_lookup_fine.mat');
% load the vss_vec. (row 21 of vss_mat)
vss_vec = vss_mat_fine(21,:);

%% initialize the SPD part
% quantization for target x values
xQuatized = -800:2:800;
% quantization for theta bin sizes
binedge = -47.5:5:47.5;
% y distance between consecutive stages
stageDist = 10;

%% initialize the MPC part
% number of steps in prediction horizon specified as prediction
% time/timeStep
nPredSteps = ceil(15/timeStep);
% bounds on heading angle
lowerBnds = -90*ones(nPredSteps,1);
upperBnds = -lowerBnds;
% initial x position
x0 = 0;
% initialization for mpc control
u0 = zeros(nPredSteps,1);
% initial tack, iniTack = {1,2}
iniTack = 1;
% dummy gate width
gate_w = 0.5;

%% the main MPC SDP algorithm
% total number of stages
numberStages = 160;
% initialize some loop parameters
minDist = 100;
% ballSize
ballSize = 0.2;
% start the SDP for loop
for ii = 0:numberStages
    if ii == 0
        xCurrent = x0;
        yCurrent = 0;
        % store position data
        posStore = [xCurrent; yCurrent];
        currentTack = iniTack;
        mpcCounter = 0;
    end
    % MPC algorithm
    % elapsed time
    elapsedTime = timeStep*mpcCounter;
    % store elapsed time
    timeStore = elapsedTime;
    % get wind at current time
    thetaCurrent = getdatasamples(windprofile,find(elapsedTime>=windprofile.time,1,'last'));
    % store wind data
    thetaStore = thetaCurrent;
    % get the correct bin for the SDP lookup table
    sdp_wind_ind = find(thetaCurrent>=binedge,1,'last');
    % find the next x location using current conditions
    xNextIdx = optimalPath(xQuatized == xCurrent,sdp_wind_ind,ii+1,iniTack);
    % next x location
    xNext = xQuatized(xNextIdx);
    % find the cost to go from the SDP cost matrix
    sdpCtg = ctg(xNextIdx,sdp_wind_ind,ii+1,:);
    % run system dynamics for 1 time step using SDP heading
    
    % run the MPC loop to reach the next point
    while minDist >= ballSize
        % print some vlaues to monitor the progression
        fprintf('postion = [%0.3f   %0.3f] m at time = %0.2f s and theta = %0.3f deg \n',...
            [xCurrent;yCurrent;elapsedTime;thetaCurrent]);
        % use PSO for coarse optimization for control sequence
        optionsPSO  = optimoptions('particleswarm','Display','off');
        [optSeqP,minFvalP] = particleswarm(...
            @(u) objfun(u,xCurrent,yCurrent,thetaCurrent,xNext,currentTack,...
            sdpCtg,vss_vec,timeStep,gate_w),nPredSteps,...
            lowerBnds,upperBnds,optionsPSO);
        % use fmincon for fine optimization for control sequence
        optionsFmincon  = optimoptions('fmincon','Display','off');
        [finSeq,finMin] = fmincon(...
            @(u) objfun(u,xCurrent,yCurrent,thetaCurrent,xNext,currentTack,...
            sdpCtg,vss_vec,timeStep,gate_w),u0,[],[],[],[],...
            lowerBnds,upperBnds,[],optionsFmincon);
        % apply the first step of the control sequence to the dynamic model
        [vel,currentTack] = boatspeed(finSeq(1),thetaCurrent,vss_vec); % total velocity and tack at every timestep
        yVel = vel*cosd(finSeq(1));    % y velocity
        yProg = yVel*timeStep;        % y progress made during each timestep
        xVel = vel*sind(finSeq(1)); % x velocity
        xProg = xVel*timeStep; % x progress
        % calculate distance from target
        minDist = norm([xNext;(ii+1)*stageDist] - [xProg;yProg]);
        % update the current location
        xCurrent = xCurrent + xProg;
        yCurrent = yCurrent + yProg;
        % increase counter
        mpcCounter = mpcCounter + 1;
        % calculate elapsed time
        elapsedTime = mpcCounter*timeStep;
        % update the measured wind
        thetaCurrent = getdatasamples(windprofile,...
            find(elapsedTime>=windprofile.time,1,'last'));
        % update data sotage matrices
        posStore = [posStore [xCurrent;yCurrent]];
        timeStore = [timeStore;elapsedTime];
        thetaStore = [thetaStore;thetaCurrent];
        
    end
    
    keyboard
    
end

%% plot scripts


