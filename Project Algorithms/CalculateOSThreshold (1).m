function [T] = CalculateOSThreshold(Pfa,N,k)
% Pfa : probability of false alarm
% N : reference window size (preferably odd if using median)
% k : number of samples in reference window of size N to use in calculating noise mean
% if only using leading or lagging window, N would represent the size of
% the leading or lagging window not the sum of the two
% outputs
% T - threshold multiplier
% PfaEst - expected Pfa for threhsold

% threshold values
possT=1:0.01:100;
if k<N/2
    k = N/2;
end
possPFA = zeros(length(possT),1);

for i=1:length(possT)
    possPFA(i) = prod((1+(possT(i)./(N+1-(1:k)))).^-1);
end
T = interp1(possPFA,possT,Pfa);


% equation for pfa as a function of T
%PfaEst = (factorial(N)/factorial(N-k))*(gamma(N-k+T+1)/gamma(N+T+1));
%semilogy(a_set,P);grid