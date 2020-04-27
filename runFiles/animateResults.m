clear
clc
close all

%% load data file
load('vidTestData.mat');
% plot update time step
samplingRate = 2;

% run animatation function
FF = animateStuff(posMPC,posSDP,windprofile,...
    'linewidth',2,...       % line widths of the boat path
    'fontsize',12,...       % font size
    'xLim',800*[-1 1],...   % x limits
    'ylim',[0 1700],...     % y limits
    'samplingRate',samplingRate);      


%% make video
% video name and format
video = VideoWriter('testVid','Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');

% video playback speed up
speedUpFactor = 50;
% set video frame rate
video.FrameRate = speedUpFactor*(1/samplingRate);

% make the video
set(gca,'nextplot','replacechildren');
open(video)
for ii = 1:length(FF)
    writeVideo(video, FF(ii));
end
close(video)