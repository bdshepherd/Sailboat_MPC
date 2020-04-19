function [optimalPath,ctg] = markov_wind_sdp_function(l_t,var)

% SDP Simple Implementation
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter', 'latex');
% set(groot, 'DefaultTextInterpreter', 'latex');
% set(groot, 'defaultAxesFontSize', 16);
% set(groot, 'defaultLineLineWidth', 1)
% set(groot,'defaultfigurecolor','w')

% clear;
% tic;
numStages = 160;   %number of stages
stageDist = 10;
yDisc = 2;
yPos = -800:yDisc:800;
yIndices = 1:length(yPos);   %choices for movement in next stage (number of positions relative to current position)
choices = yIndices;
wangChoices = -45:5:45;    %possible wind angles
wangIndices = 1:length(wangChoices);
tackPenalty = 5; %sec

lengthScale_pred = l_t; %min
overallVariance_pred = var;

pathDist = sqrt(stageDist^2 + (yPos - yPos').^2);
boatAngle = -atand((yPos - yPos')/stageDist);  %angle to go from y position 1 to y position 2
tmpWang(1,1,:) = wangChoices;
twa = boatAngle - tmpWang;
%% velpol here
load('vss_lookup_fine.mat'); % rows: 0:.1:10 m/s windspeed, columns: 0:.1:180 TWA 
                             % Use 2 m/s wind for all simulations (row 21)
                             % Note 0 m/s boat velocities modified to very
                             % low values to avoid infinite ctg
                             % calculations
windSpeed = 2;
dWindSpeed = .1;
spdIndex = 1 + floor(windSpeed/dWindSpeed);
% angIndex = 1 + ceil(abs(apparentWindAngle)*(length(vss_mat_fine)-1)/180);
angIndex = 1 + round(abs(twa)*10);
boatSpeed = reshape(vss_mat_fine(spdIndex,angIndex),size(twa,1),size(twa,2),size(twa,3));
%% rest of code
stageCostMat = pathDist./boatSpeed;

ctg = zeros(length(yPos),length(wangChoices),numStages,2);
optimalPath = ctg;

%% determining markov matrix

timestep = (13/80)*stageDist/6; %approx. timestep in minutes based on grid size
K_pred = timestep/lengthScale_pred;

w = sqrt(overallVariance_pred*K_pred/(.42 + K_pred));

markTrans_pred = markovSDP(K_pred,w^2,wangChoices);

%% calculating cost to go/optimal path using SDP
% DO NOT NEED TO RUN THIS! 
% Output for various lt and overall variance found in the
% 'SDP_lookup_tables' folder'


% count = 0;
% update = 1;
% for i = numStages:-1:1  %for each stage
%     if mod(i,update) == 0
%         count = count+1;
%         disp(i);
%         if count~=1
%             timestg = seconds(toc);
%             esttimerem = i*timestg/update;
%             disp(['Expected remaining time: ',datestr(esttimerem,'HH:MM:SS')])
%             tic
%         else
%             tic;
%         end
%     end
%     if i == numStages; choices = floor(length(yPos)/2-3):floor(length(yPos)/2+3); else; choices = yIndices; end
%     for j = yIndices  %for each y position
%         for k = wangIndices   %for each possible wind angle
%             for m = choices   %for each path
%                 for n = [1 2]
%                     
%                     %                 wang = wangChoices(k);
%                     %                 boatSpeed = abs(sind(boatAngle(j,m) - wang)) + .0001;
%                     %                 stageCost = pathDist(j,m)/boatSpeed;  %time to reach next node
%                     %absolute value of sine is estimate of velocity polar
%                     nextAng = atand((m-j)*yDisc/stageDist);
%                     nextTack = 1 + (nextAng > wangChoices(k));
%                     if n == nextTack
%                         q = 0;
%                     else
%                         q = tackPenalty;
%                     end
%                     
%                     stageCost = stageCostMat(j,m,k) + q;
%                     
%                     %calculating cost to go instead of stage cost
%                     if i ~= numStages
%                         windDist = (wangIndices==k)*markTrans_pred;  %wind distribution
%                         cost = stageCost + windDist*ctg(m,:,i+1,nextTack)';  %including cost to go at (expected) next node
%                     else
%                         cost = stageCost;   %if it is the final node, cost to go is just the stage cost
%                     end
%                     
%                     if ctg(j,k,i,n) == 0 || cost < ctg(j,k,i,n)   %store cost to go if it is the best one
%                         ctg(j,k,i,n) = cost;
%                         optimalPath(j,k,i,n) = m;
%                     end
%                 end
%             end
%         end
%     end
% end


%% Plotting optimal path from set start position

% % Note: for determining wind only, you can use something similar to the 
% %   following. You would need to keep track up the wind angle. rng(i) is
% %   used to keep the runs consistent so they can be recreated. First wind
% %   angle is always zero
% %   
% for i = 1:numPlots
%     wind = (length(wangChoices)+1)/2 ;
%     rng(i);
%     for j = 1:numStages
%        weights = markTrans_act(wind,:);   
%        wind = randsample(1:length(wangChoices),1,true,weights);         
%     end    
% end



% Missing a few lines here for selecting the volatility parameters you are
% interested in

tf = 0;
if tf
filename = ['sdp_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act),'.mat'];
load(filename);
ctg = eval(['ctg_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act)]);
optimalPath = eval(['optpath_lt',num2str(lengthScale_act),'_ovar',num2str(overallVariance_act)]);

timestep = ctg(ceil(length(yPos)/2),(length(wangChoices)+1)/2,1,1)/(numStages*60); 
K_act = timestep/lengthScale_act;
w_act = sqrt(overallVariance_act*K_act/(.42 + K_act));
markTrans_act = markovSDP(K_act,w_act^2,wangChoices);

numPlots = 2000;
figure;
hold on
axis([0 numStages*stageDist yPos(1) yPos(end)])
times = zeros(1,numPlots);
path = zeros(numPlots,numStages+1);
sdp_numtacks = zeros(numPlots,1);
for i = 1:numPlots
    wind = (length(wangChoices)+1)/2 ;
    path(i,1) = ceil(length(yPos)/2);
    time = zeros(1,numStages);
    currTack = 1;
    rng(i);
    for j = 1:numStages
        path(i,j+1) = optimalPath(path(i,j),wind,j,currTack);
        time(j) = stageCostMat(path(i,j),path(i,j+1),wind);
        boatDir = atand((path(i,j+1)-path(i,j))*yDisc/stageDist);
        newTack = 1 + (boatDir > wangChoices(wind));
        if newTack ~= currTack
            time(j) = time(j) + tackPenalty;
            sdp_numtacks(i) = sdp_numtacks(i) + 1;
        end
        currTack = newTack;
        weights = markTrans_act(wind,:);        % u
        wind = randsample(1:length(wangChoices),1,true,weights);
    end
    tmp = plot(stageDist*(0:numStages),yPos(path(i,:)),'r-',stageDist*(0:numStages),-yPos(path(i,:)),'r-');
    tmp(1).Color(4) = .2;
    tmp(2).Color(4) = .2;
    times(i) = sum(time);
%     if max(time)>200
%         figure
%         plot(time)
%         error('Long time')
%     end
end

percentiles = sort([yPos(path); -yPos(path)]);
bound90 = percentiles(floor(.9*numPlots*2),:);
bound10 = percentiles(ceil(.1*numPlots*2),:);
bound30 = percentiles(ceil(0.3*numPlots*2),:);
bound70 = percentiles(floor(.7*numPlots*2),:);
bnds10 = plot(stageDist*(0:numStages),bound10,'k','DisplayName','$P_{10/90}$');
bnds90 = plot(stageDist*(0:numStages),bound90,'k');
bnds30 = plot(stageDist*(0:numStages),bound30,'--k','DisplayName','$P_{30/70}$');
bnds70 = plot(stageDist*(0:numStages),bound70,'--k');
% bnds30 = plot(stageDist*(0:numStages),bound30,'--k',stageDist*(0:numStages),bound70,'--k');
bnds10.LineWidth = 2;
bnds30.LineWidth = 2;
bnds70.LineWidth = 2;
bnds90.LineWidth = 2;
legend([bnds10 bnds30],'Location','Southeast')
title(['SDP results when $l_t=',num2str(lengthScale_act),'$ and $\sigma^2=',num2str(overallVariance_act),'$'])
ylabel('X Position (m)')
xlabel('Y Position (m)')
set(gca,'View',[90 -90])
% text(50,-400,['SDP Expected Time = ',num2str(round(ctg(path(1),(length(wangChoices)+1)/2,1),2))])
% text(50,-450,['Simulated Average Time = ',num2str(round(mean(times),2))])
hold off
title(['SDP results when lt =',num2str(lengthScale_act),' and ovar = ',num2str(overallVariance_act)])
% figure
% plot(times)

% figure
% histogram(times-ctg(path(1),(length(wangChoices)+1)/2,1))
% xlabel('Difference between actual and expected cost to go (sec)')
% title(['Predicted $l_t=',num2str(lengthScale_pred),'$ and $\sigma^2=',num2str(overallVariance_pred),'$ vs actual $l_t=',num2str(lengthScale_act),'$ and $\sigma^2=',num2str(overallVariance_act),'$'])

%fprintf('Predict length scale = %f & overall variance = %f, but actual length scale = %f with overall variance %f\n',lengthScale_pred,overallVariance_pred,lengthScale_act,overallVariance_act)
fprintf('Simulated Average Time: %f\n', mean(times))
fprintf('SDP Expected Cost to Go: %f\n', ctg(path(1),(length(wangChoices)+1)/2,1))
%toc 
fprintf('\n')

%     end
% end
else
    timestep = exp_time/(numStages*60); %approx. timestep in minutes (should be updated based on velpol)
    K_act = timestep/lengthScale_act;
    w_act = sqrt(overallVariance_act*K_act/(.42 + K_act));
    markTrans_act = markovSDP(K_act,w_act^2,wangChoices);
end

%%
% FolderName = 'C:\Users\bshephe\Google Drive\Graduate Research\MATLAB and Simulink\pred_lt_10_ovvar_3500';   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   savefig(FigHandle, fullfile(FolderName, [FigName, '.fig']));
% end
% 



%440 sec for 121 yPos, 300 stages, 15 wind directions
% total is 1237 calculations/sec (some efficiency gains still planned)
%savefig(['LenSc',num2str(lengthScale),'_OvrVar',num2str(overallVariance),'_Coarse.fig'])

end