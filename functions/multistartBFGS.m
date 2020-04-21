function [optDsgn,minF] = multistartBFGS(objF,X,lb,ub,varargin)
%PARTICLESWARMOPT
% Function to find the global minimum of an objective function using
% multi start BFGS Optimization, which is basically particle swarm
% optimization with a smarter optimum seeking
%
% Required inputs:
% objF = objective function which should be a function of X
% X =  dummy design matrix, rows represent design variables and columns represent different
% designs. Can also be the initial condition
% lb, ub =  lower and upper bounds on the design, should be the same size
% as X
%
% varargin parameters:
% 'nStarts' = number of starting points

% Example use of code
% [maxF,optDsgn] = multistartBFGS(@(x)objF(x),X,lb,ub,'nStarts',20);

p = inputParser;
addRequired(p,'objF');
addRequired(p,'X',@isnumeric);
addRequired(p,'lb',@isnumeric);
addRequired(p,'ub',@isnumeric);
addParameter(p,'nStarts',25,@isnumeric);

parse(p,objF,X,lb,ub,varargin{:});

ss = p.Results.nStarts;

% design space size
dsgnSize = size(X);
% initial swarm size
iniSwarm = lb + (ub-lb).*rand([dsgnSize,ss]);
swarm = NaN(size(iniSwarm));

% initial design space
fVal = NaN(ss,1);
options  = optimoptions('fmincon','Display','off');

for ii = 1:ss
    
    [OptPts,OptFval] = fmincon(@(X) p.Results.objF(X),iniSwarm(:,:,ii),...
        [],[],[],[],lb,ub,[],options);
    
    fVal(ii) = OptFval;
    swarm(:,:,ii) = OptPts;
    
end

[minF,idx] = min(fVal);
optDsgn = swarm(:,:,idx);

end





