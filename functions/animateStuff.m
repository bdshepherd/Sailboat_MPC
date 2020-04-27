function aa = animateStuff(mpcPosTs,sdpPosTs,varargin)

%% parse the input
p = inputParser;
% add a required time series object for positions using MPC
addRequired(p,'mpcPos', @(x) isa(x,'timeseries'));
% add a required time series object for positions using SDP
addRequired(p,'sdpPos', @(x) isa(x,'timeseries'));
% add optional line width parameter
addParameter(p,'linewidth',1,@isnumeric);
% add optional font size parameter
addParameter(p,'fontSize',12,@isnumeric);
% add optimal x limits
addParameter(p,'xLim',800*[-1 1],@isnumeric);
% add optimal y limits
addParameter(p,'yLim',1700*[0 1],@isnumeric);
% add optimal video file name
addParameter(p,'vidFileName','testVid',@ischar);
% add optimal video file name
addParameter(p,'frameRate',20,@isnumeric);
% data plotting rate
addParameter(p,'samplingRate',2,@isnumeric);

% parse the inputs
parse(p,mpcPosTs,sdpPosTs,varargin{:});

%% data processing
destinationIdx = find(abs(mpcPosTs.Data(:,2)-1600)<0.1,1);

% resample data
tVec = 0:p.Results.samplingRate:p.Results.mpcPos.Time(destinationIdx);
mpcPosTs = resample(p.Results.mpcPos,tVec);
% number of time steps
nSteps = numel(tVec);
% extract position vectors from mpc
mpcPosVals = squeeze(mpcPosTs.Data);
% resample sdp at the same time vector
sdpPosTs = resample(p.Results.sdpPos,tVec);
% extract position vectors from sdp
sdpPosVals = squeeze(sdpPosTs.Data);
% rearrage MPC the data matrices
mpcReIdx = size(mpcPosVals) == nSteps;
mpcPosVals = permute(mpcPosVals,cat(2,find(mpcReIdx),find(~mpcReIdx)));
% rearrage SDP the data matrices
sdpReIdx = size(sdpPosVals) == nSteps;
sdpPosVals = permute(sdpPosVals,cat(2,find(sdpReIdx),find(~sdpReIdx)));


%% plot the data
% custome color scheme
colorVals = colorScheme;
% open a figure
figure;
set(gca,'FontSize',p.Results.fontSize);
% preallocate the get frame structure
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);
% start the for loop for each tVal
for ii = 1:nSteps
    
    % set axis properties at the first iteration and delete things for the
    % rest
    if ii == 1
        hold on
        grid on
        xlim(p.Results.xLim);
        ylim(p.Results.yLim);
        xlabel('x (m)');
        ylabel('y (m)');
    else
        h = findall(gca,'type','line');
        delete(h);
    end
    
    % plot the MPC positions
    plot(mpcPosVals(1:ii,1),mpcPosVals(1:ii,2),...
        'linewidth',p.Results.linewidth,...
        'color',colorVals(1,:));
    % plot the SDP positions
    plot(sdpPosVals(1:ii,1),sdpPosVals(1:ii,2),...
        'linewidth',p.Results.linewidth,...
        'color',colorVals(2,:));
    
    legend('MPC','Pure SDP','Location','northeast');
    title(sprintf('Time = %0.2f s',tVec(ii)));
    
    % get the frames
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

%% make the video
% % % % video setting
video = VideoWriter(p.Results.vidFileName,'Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = p.Results.frameRate;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)



% % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % %
    function c = colorScheme()
        c = repmat(1/255.*[228,26,28
            55,126,184
            77,175,74
            152,78,163
            255,127,0
            255,255,51],10,1);
    end





end