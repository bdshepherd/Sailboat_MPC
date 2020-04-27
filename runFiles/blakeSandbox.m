clear
clc
% Blake Sandbox for testing

%% generate stochastic wind
% length scale of stochastic wind
lengthScale_act = 10;
% variance of stochastic wind
overallVariance_act = 2000;
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
%Augment ctg with zeros for ctg from end of course to end
ctg(:,:,161,:) = 0;
% evaluate the optimal path matrix
optimalPath = eval(['optpath_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act)]);
% load the vss_mat matrix
load('vss_lookup_fine.mat');
% load the vss_vec. (row 21 of vss_mat)
vss_vec = vss_mat_fine(21,:);
%For generating 1-tack solutions that may work


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
for ii = 1:numberStages
    if ii == 1 %first stage
        xCurrent = x0;
        yCurrent = 0;
        % store position data
        posStore = [xCurrent yCurrent;zeros(36000,2)];
        currentTack = [iniTack;zeros(36000,1)];
        mpcCounter = 0;
        timeCurrent = 0;
        timeStore = zeros(36001,1); % get wind at current time
        thetaCurrent = getdatasamples(windprofile,find(timeCurrent>=windprofile.time,1,'last'));
        % store wind data
        thetaStore = [thetaCurrent;zeros(36000,1)];
        % store time of next wind update
        timeupdate.wind = windprofile.time(find(windprofile.time>timeCurrent,1));
        tackStore = zeros(36001,1);
        uStore = zeros(36001,1);
        sdp_last_update = 0;
    end
    % MPC algorithm
%     % elapsed time
%     elapsedTime = timeStep*mpcCounter;
%     % store elapsed time
%     timeStore = elapsedTime;
   sdp_last_update = yCurrent;
    % get the correct bin for the SDP lookup table
    sdp_wind_ind = find(thetaCurrent>=binedge,1,'last');
    % find nearest x position
    [~,xind] = min(abs(xCurrent - xQuantized));
    % find the next x location using current conditions
    xNextIdx = optimalPath(xind,sdp_wind_ind,ii,iniTack);
    % next x location
    xNext = xQuantized(xNextIdx);
    % find the cost to go from the SDP cost matrix
    sdpCtg = ctg(xNextIdx,sdp_wind_ind,ii+1,:);
    % run system dynamics for 1 time step using SDP heading
    
    % run the MPC loop to reach the next point
    cond2 = true;
    sdp_trigger = false;
    while ~sdp_trigger && cond2  % get 98% to next SDP stage
        % Reset sdp_trigger
        sdp_trigger = false;
        % Change horizon length based on y-distance from next waypoint
        nPredSteps = 0+ceil(4*(ii*stageDist - yCurrent))/timeStep; % 25 seconds of prediction
        % bounds on heading angle
        lowerBnds = -90*ones(nPredSteps,1);
        upperBnds = -lowerBnds;
        heading2waypoint = atan2d(xNext - xCurrent,ii*stageDist - yCurrent);
%%
        u_best = zeros(1,nPredSteps); %init
        f_best = 999999; % init        
        psi_all = 37:.05:39; % Based on VMG analysis
        for j = psi_all
            % Based on general VMG
            psi_max = thetaCurrent + j;
            psi_min = thetaCurrent - j;
            % temp matrices for headings
            port_mat = psi_max*ones(nPredSteps);
            star_mat = psi_min*ones(nPredSteps);
            % Create matrix of potential initial paths
            init_guess_mat = [heading2waypoint*ones(1,nPredSteps);...
                tril(port_mat) + triu(star_mat,1); tril(star_mat) + triu(port_mat,1)];
            % Find direct or one-tack solution based on heading angle limits
            % with smallest objfun value
            
            for i = 1:size(init_guess_mat,1) % Loop over all initial guesses
                objF = objfun(init_guess_mat(i,:),xCurrent,yCurrent,thetaCurrent,xNext,ii*stageDist,currentTack(1),...
                    sdpCtg,vss_vec,timeStep,gate_w);
                if objF<f_best
                    f_best = objF; %replace best objfun value
                    u_best = init_guess_mat(i,:); % replace best found u
                end
            end
        end
%%
%         % use PSO for coarse optimization for control sequence
%         optionsPSO  = optimoptions('particleswarm',...
%             'SwarmSize',200,'UseParallel',true,'MaxStallIterations',10,...
%             'InitialSwarmSpan',180,'InitialSwarmMatrix',heading2waypoint*ones(1,nPredSteps),'Display','off');
%         [optSeqP,minFvalP] = particleswarm(...
%             @(u) objfun(u,xCurrent,yCurrent,thetaCurrent,xNext,ii*stageDist,currentTack(1),...
%             sdpCtg,vss_vec,timeStep,gate_w),nPredSteps,...
%             lowerBnds,upperBnds,optionsPSO);
%         optSeqP
%         use fmincon by itself
%%
% %Fast option that selects the best of all the possible options tested
% above
        finSeq = u_best;
        
% %         Fmincon optimization 
%         optionsFmincon  = optimoptions('fmincon','UseParallel',true,'MaxFunctionEvaluations',1000,'Display','off');
%         [finSeq,finMin] = fmincon(...
%             @(u) objfun(u,xCurrent,yCurrent,thetaCurrent,xNext,ii*stageDist,currentTack(1),...
%             sdpCtg,vss_vec,timeStep,gate_w),u_best,[],[],[],[],...
%             lowerBnds,upperBnds,[],optionsFmincon);

        % store time of next MPC update
        timeupdate.mpc = timeCurrent + timeStep;
        % dynamic modelling     
        cond3 = true;
        while timeCurrent<timeupdate.mpc && cond3 %not yet time for MPC update
            % print some vlaues to monitor the progression
            fprintf('stage = %d,pos = [%0.3f   %0.3f] m, time = %0.2f s,theta = %0.3f deg, u = %0.3f deg \n',...
                [ii;xCurrent;yCurrent;timeCurrent;thetaCurrent;finSeq(1)]);
            % Under assmuption that wind and control sequence won't change
            [vel,currentTack] = boatspeed(finSeq,thetaCurrent,vss_vec); % total velocity and tack at every timestep
            tackStore(storeIdx) = currentTack(1); % store tack
            uStore(storeIdx) = finSeq(1); % store control decision
            yVel = vel.*cosd(finSeq);    % y velocity
            yProg = yVel*timeStep; % y progress made during each timestep
            yPos = yCurrent + cumsum(yProg); % position after each timestep
            xVel = vel.*sind(finSeq); % x velocity
            xProg = xVel*timeStep; % x progress
            xPos = xCurrent + cumsum(xProg); % position after each timestep
            if yPos(1)>=ii*stageDist %Can reach next SDP stage before next MPC update under current wind
                ts_end_per = (ii*stageDist - yCurrent)/yProg(1); %proportion of MPC timestep needed
                timeupdate.sdp = timeCurrent + ts_end_per*timeStep;
            else % Doesn't reach end of sdp stage before mpc update
                timeupdate.sdp = timeCurrent + 100; % Set far in future
            end
            
            %Find appropriate deltaT for dynamics
            [timeupdate.next,timeupdate.index] = min([timeupdate.sdp timeupdate.mpc timeupdate.wind]);
            switch timeupdate.index
                case 1 %sdp update would happen first
                    deltat = ts_end_per*timeStep;
                    sdp_trigger = true;
                case 2 % MPC update would happen first
                    deltat = timeupdate.mpc - timeCurrent;
                case 3 % wind update comes first
                    deltat = timeupdate.wind - timeCurrent;
            end
            % Use portion of dynamics already calculated above
            xCurrent = xCurrent + xVel(1)*deltat; % Position update
            yCurrent = yCurrent + yVel(1)*deltat;  % Position update
            timeCurrent = timeCurrent + deltat;  % Time update
            % update the measured wind
            thetaCurrent = getdatasamples(windprofile,...
                find(timeCurrent>=windprofile.time,1,'last'));
            % store time of next wind update
            timeupdate.wind = windprofile.time(find(windprofile.time>timeCurrent,1));
            % update data sotage matrices
            storeIdx = storeIdx + 1; % update index for storage
            timeStore(storeIdx) = timeCurrent;
            posStore(storeIdx,1:2) = [xCurrent yCurrent];           
            thetaStore(storeIdx) = thetaCurrent;
            if ii==numberStages
                cond2 = true;
            elseif (yCurrent - sdp_last_update)<8 %prevent back-to-back updates
                cond2 = true;
                sdp_trigger = false;
            else
                cond2 = mod(yCurrent,stageDist)/stageDist <=.99;
            end
            cond3 = ~sdp_trigger && cond2;
        end
    end
    
%     keyboard
    
end

%%
% Calculate total tacking penalties
tackPen = 5; % [sec] per each tack across wind
numtacks = nnz(diff(tackStore)); % total number of tacks across wind
timeTackTot = tackPen*numtacks;
end_sim_idx = find(timeStore,1,'last');
timeStore = [timeStore(1:end_sim_idx-1); timeStore(end_sim_idx)+timeTackTot];
% Truncate stored variables and create timeseries
posStore = posStore(1:end_sim_idx,:);
uStore = uStore(1:end_sim_idx);
thetaStore = thetaStore(1:end_sim_idx);
% Create timeseries' (see collection at bottom)
pos = timeseries(posStore,timeStore,'Name','Position');
pos.DataInfo.Units = 'meters';
pos.TimeInfo.Units = 'seconds';

useq = timeseries(uStore,timeStore,'Name','Heading Control Sequence');
useq.DataInfo.Units = 'degrees';
useq.TimeInfo.Units = 'seconds';
useq.DataInfo.Interpolation = 'zoh';

theta = timeseries(thetaStore,timeStore,'Name','Theta');
theta.DataInfo.Units = 'degrees';
theta.TimeInfo.Units = 'seconds';
theta.DataInfo.Interpolation = 'zoh';
% Timeseries Collection
tsc = tscollection({pos useq theta},'Name','Run Data');










