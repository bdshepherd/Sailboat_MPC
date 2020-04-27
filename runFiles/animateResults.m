clear
clc

%% load data file
load('vidTestData.mat');

% run animatation function
FF = animateStuff(posMPC,posSDP,windprofile);

%% make video
video = VideoWriter('testVid','Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 20;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(FF)
    writeVideo(video, FF(ii));
end
close(video)