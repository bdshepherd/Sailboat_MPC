function P = markovSDP(K,w,psi_all)
%markovSDP generates a Markov transition probability matrix for running
%stochastic dynamic programming simulations according to 
% psi(k+1) = psi(k) + K*(psi_prevailing - psi(k)) + w(k)
%   Inputs:
%           K = number between 0 and 1 to represent the jump scaling
%           between current state and the mean value of the pdf for the
%           next state
%           w = variance of normal distributed random variable
%           psi_all = vector of possible wind directions
%           (should be odd for symmetry and have equal bin sizes)
% Outputs:
%           P = probability distribution matrix 

if mod(length(psi_all),2)==0
    error('Number of states in the Markov model must be odd for symmetry')
end

binw = psi_all(2)-psi_all(1); %bin width
binedge = psi_all(1)-binw/2:binw:psi_all(end)+binw/2; % Bin Edges
P = zeros(length(psi_all)); % Initialize P

for i = 1:length(psi_all)
    pd = makedist('Normal','mu',psi_all(i)*(1-K),'sigma',w^(1/2)); %Make new distribution shifting mean and using input variance
    edgecdf = cdf(pd,binedge); % CDF and bin edges
    P(i,:) = edgecdf(2:end)-edgecdf(1:end-1); %PDF for individual bins
    if sum(P(i,:))<1 %If row does not have a total probability of 1
         P(i,:) = P(i,:)./sum(P(i,:)); % Normalize probabilities to have sum of 1
    end    
end

end

