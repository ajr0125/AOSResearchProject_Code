function [pd] = CMLD_CFAR(numTrials, PFA, N, magSqGaussNoise, signalData, numSampsCensored)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% T - threshold multiplier
% SNR_dB - vector of SNR values in dB

% pd - probability of detection
numHomogeneousSamples = N - numSampsCensored;
sortedmagSqGaussNoise = sort(magSqGaussNoise);
T = (PFA.^(-1/(numHomogeneousSamples-1)))-1;
cmldNME =  (sortedmagSqGaussNoise(numHomogeneousSamples,:) * (N-numHomogeneousSamples) + ...
    sum(sortedmagSqGaussNoise(1:numHomogeneousSamples-1,:)));
pd = sum(signalData > T .* cmldNME)/numTrials;
