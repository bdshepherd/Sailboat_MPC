function pl = plotStuff(tsObject,varargin)

%% parse the input
p = inputParser;
% add a required time series object
addRequired(p,'tsObject', @(x) isa(x,'timeseries'));
% add optional line width parameter
addParameter(p,'linewidth',1,@isnumeric);
% add optional font size parameter
addParameter(p,'fontSize',12,@isnumeric);
% add optional font size parameter
addParameter(p,'ylabel','',@ischar);

% parse the inputs
parse(p,tsObject,varargin{:});

% extract time data
tsTime = squeeze(p.Results.tsObject.Time);
% extract data vals
tsData = squeeze(p.Results.tsObject.Data);
% rearrage the data matrix
szData = size(tsData);
reIdx = szData == numel(tsTime);
tsData = permute(tsData,cat(2,find(reIdx),find(~reIdx)));

%% plot
% custome color scheme
colorVals = colorScheme;
% plot the data
for ii = 1:szData(szData~=numel(tsTime))
    pl = plot(tsTime,tsData(:,ii),...
        'linewidth',p.Results.linewidth,...
        'color',colorVals(ii,:));
    hold on
end
grid on
xlabel('Time (sec)')
ylabel(p.Results.ylabel);

set(findobj('-property','FontSize'),'FontSize',p.Results.fontSize)


    function c = colorScheme()
        c = 1/255.*[228,26,28
            55,126,184
            77,175,74
            152,78,163
            255,127,0
            255,255,51];
    end


end


