clear
clc
% Blake MPC investigation

%% generate stochastic wind
% length scale of stochastic wind
lengthScale_act = 120;
% variance of stochastic wind
overallVariance_act = 50;
% time step for wind
tsw = 0.2;
% final time for wind gereration as wind it is generated as a timeseries
tFinal = 60*60; % minutes*60
% random number generator seed for stochastic wind generation
rngSeed = 1;
% run the stochastic wind generation function
windprofile = windprofileGen(tsw,lengthScale_act,overallVariance_act,...
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
xQuantized = -800:2:800;
% quantization for theta bin sizes
binedge = -47.5:5:47.5;
% y distance between consecutive stages
stageDist = 10;

%% initialize the MPC part
% time step for MPC
timeStep = 1;
% initial x position
x0 = 0;
% initialization for mpc control
% u0 = zeros(nPredSteps,1);
% initial tack, iniTack = {1,2}
iniTack = 2;
% dummy gate width
gate_w = .2;

%% the main MPC SDP algorithm
% total number of stages
numberStages = 160;
% initialize some loop parameters
% minDist = 100;
% % ballSize
% ballSize = 0.2;
% Initialize variable for indexing variable storage
storeIdx = 1;
% start the SDP for loop

%% Deviation from run script to investigate wind from different angles
% Note that the waypoint is fixed as if you were in a direct headwind, but
% the wind angle is actually the outer for loop
ii = 1;
thetaVals = 0:3:45;

currentTack = 2;
xCurrent = x0;
yCurrent = 0;
timeCurrent = 0;
% find the next x location using set conditions
xNextIdx = optimalPath(401,10,1,2);
% next x location
xNext = xQuantized(xNextIdx);
% find the cost to go from the SDP cost matrix
sdpCtg = ctg(xNextIdx,10,ii+1,:);
% optimization variables
% Change MPC prediction horizon length
nPredSteps = 0+ceil(3*(ii*stageDist - yCurrent))/timeStep;
% bounds on heading angle
lowerBnds = -90*ones(nPredSteps,1);
upperBnds = -lowerBnds;

% provide direct path there control sequence
heading2waypoint = atan2d(xNext - xCurrent,ii*stageDist - yCurrent);

%For generating 1-tack solutions that may work
psi_all = 38; % Based on VMG analysis

% Initialize storage variables
    % rows corresponds to thetaVals above
    % ps is for particle swarm
    % fmin is if particle swarm with fmincon refinement

ps_u = zeros(length(thetaVals),nPredSteps); % ps control sequences
ps_x = zeros(length(thetaVals),nPredSteps+1); % ps x positions
ps_y = zeros(size(ps_x)); % ps y positions

fmin_u = zeros(length(thetaVals),nPredSteps); % fmin control sequences
fmin_x = zeros(length(thetaVals),nPredSteps+1); % fmin x positions
fmin_y = zeros(size(fmin_x)); % fmin y positions
tic
for j = 1:length(thetaVals)
 
    thetaCurrent = thetaVals(j);
    fprintf('Working on wind angle theta = %0.2f \n',thetaCurrent)
    % Based on general VMG
    psi_max = thetaCurrent + psi_all;
    psi_min = thetaCurrent - psi_all;
    % temp matrices for headings
    port_mat = psi_max*ones(nPredSteps);
    star_mat = psi_min*ones(nPredSteps);
    % Create matrix of potential initial paths
    init_guess_mat = [heading2waypoint*ones(1,nPredSteps);...
        tril(port_mat) + triu(star_mat,1); tril(star_mat) + triu(port_mat,1)];
    

    % use PSO for coarse optimization for control sequence
%     optionsPSO  = optimoptions('particleswarm',...
%         'SwarmSize',200,'UseParallel',true,'MaxStallIterations',10,...
%         'InitialSwarmSpan',180,'InitialSwarmMatrix',init_guess_mat,'Display','off');
%     [optSeqP,minFvalP] = particleswarm(...
%         @(u) objfun(u,xCurrent,yCurrent,thetaCurrent,xNext,ii*stageDist,currentTack(1),...
%         sdpCtg,vss_vec,timeStep,gate_w),nPredSteps,...
%         lowerBnds,upperBnds,optionsPSO);
    
    % Fmincon alone temporarily
    u_best = zeros(1,nPredSteps);
    f_best = 999999;
    for i = 1:size(init_guess_mat,1)
        objF = objfun(init_guess_mat(i,:),xCurrent,yCurrent,thetaCurrent,xNext,ii*stageDist,currentTack(1),...
            sdpCtg,vss_vec,timeStep,gate_w);
        if objF<f_best
            f_best = objF;
            u_best = init_guess_mat(i,:);
        end
    end
        %         use fmincon for fine optimization for control sequence
        optionsFmincon  = optimoptions('fmincon','UseParallel',true,'MaxFunctionEvaluations',3000,'Display','off');
        [finSeq,finMin] = fmincon(...
            @(u) objfun(u,xCurrent,yCurrent,thetaCurrent,xNext,ii*stageDist,currentTack(1),...
            sdpCtg,vss_vec,timeStep,gate_w),u_best,[],[],[],[],...
            lowerBnds,upperBnds,[],optionsFmincon);
      
%     % Using only particle swarm solution
%     [vel,currentTack] = boatspeed(optSeqP,thetaCurrent,vss_vec); % total velocity and tack at every timestep
%     ps_u(j,:) = optSeqP;
%     yVel = vel.*cosd(optSeqP);    % y velocity
%     yProg = yVel*timeStep; % y progress made during each timestep
%     ps_y(j,:) = [yCurrent [yCurrent + cumsum(yProg)]]; % position after each timestep
%     xVel = vel.*sind(optSeqP); % x velocity
%     xProg = xVel*timeStep; % x progress
%     ps_x(j,:) = [xCurrent [xCurrent + cumsum(xProg)]]; % position after each timestep 
    
    % Particle Swarm with Fmincon refinement path
    [vel,currentTack] = boatspeed(finSeq,thetaCurrent,vss_vec); % total velocity and tack at every timestep
    tackStore(storeIdx) = currentTack(1); % store tack
    uStore(storeIdx) = finSeq(1); % store control decision
    fmin_u(j,:) = finSeq;
    yVel = vel.*cosd(finSeq);    % y velocity
    yProg = yVel*timeStep; % y progress made during each timestep
    fmin_y(j,:) = [yCurrent [yCurrent + cumsum(yProg)]]; % position after each timestep
    xVel = vel.*sind(finSeq); % x velocity
    xProg = xVel*timeStep; % x progress
    fmin_x(j,:) = [xCurrent [xCurrent + cumsum(xProg)]];  % position after each timestep    
    
end
toc

    
%% plot scripts

% Sample plotting
figure;
hold on
for i = 1:size(ps_u,1)
    tmp = plot(ps_x(i,:),ps_y(i,:),'r-');
    tmp.Color(4) = .2;
    
end
plot(xNext,ii*stageDist,'*k','MarkerSize',8)
title('Particle Swarm MPC paths For all Wind Angles')

figure;
hold on
for i = 1:size(ps_u,1)
    tmp2 = plot(fmin_x(i,:),fmin_y(i,:),'b--');
    tmp2.Color(4) = .2;    
end
plot(xNext,ii*stageDist,'*k','MarkerSize',8)
title('Particle Swarm wint Fmincon Refinement MPC paths For all Wind Angles')

