clear
clc
% run file to test the various functions

%% initialization
% load the vss_mat matrix
load('vss_lookup_fine.mat');
% load the vss_vec. (row 21 of vss_mat)
vss_vec = vss_mat_fine(21,:);
% dummy starting x location
x = 0;
% dummy starting y location
y = 1; 
% dummy starting theta value
theta = 0; 
% dummy target x location
xd = 2; 
% startin tack
k = 1; 
% length of prediction horizon
nPred = 50; 
% control sequence. Vector of number of stages
u = repmat(atan((xd-x)/10)*180/pi,nPred,1); 
% dummy gate width
gate_w = 0.5;
% dummy time step
ts = 0.1;

%% find the respective cost to go values
% length scale
lengthScale_pred = 10; 
% prediction variance
overallVariance_pred = 500;
% create file name using length scale and prediction variance
filename = ['sdp_lt',num2str(lengthScale_pred),'_ovar',num2str(overallVariance_pred),'.mat'];
% load the respective file
load(filename);
% evaluate the cost to go matrix
ctg = eval(['ctg_lt',num2str(lengthScale_pred),'_ovar',num2str(overallVariance_pred)]);
% evaluate the optimal path
optimalPath = eval(['optpath_lt',num2str(lengthScale_pred),'_ovar',num2str(overallVariance_pred)]);
% find the appropriate cost to go
sdp_ctg = ctg(xQuatized == xd,thetaQuantized == theta,k,:);
% test if the objective function works
J = objfun(u,x,y,theta,xd,k,sdp_ctg,vss_vec,ts,gate_w);

%% test optimization
% lower bounds
lowerBnds = -45*ones(nPred,1);
upperBnds = 45*ones(nPred,1);

[optSeqPSO,minFvalPSO,allData] = particleSwarmMinimization(...
    @(u) objfun(u,x,y,theta,xd,k,sdp_ctg,vss_vec,ts,gate_w),u,...
    lowerBnds,upperBnds,'swarmSize',40,'cognitiveLR',0.4,...
    'socialLR',0.2,'maxIter',100);

[optSeqMBFGS,minFvalMBFGS] = multistartBFGS(...
    @(u) objfun(u,x,y,theta,xd,k,sdp_ctg,vss_vec,ts,gate_w),u,...
    lowerBnds,upperBnds,'nStarts',20);

[optSeqP,minFvalP] = particleswarm(...
    @(u) objfun(u,x,y,theta,xd,k,sdp_ctg,vss_vec,ts,gate_w),nPred,...
    lowerBnds,upperBnds);

options  = optimoptions('fmincon','Display','off');
[fimSeq,finMin] = fmincon(...
    @(u) objfun(u,x,y,theta,xd,k,sdp_ctg,vss_vec,ts,gate_w),u,[],[],[],[],...
    lowerBnds,upperBnds,[],options);



