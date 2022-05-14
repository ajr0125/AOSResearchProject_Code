function [pd] = CA_CFAR(numTrials, PFA, N, magSqGaussNoise, signalData)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% T - threshold multiplier
% SNR_dB - vector of SNR values in dB

% pd - probability of detection
T = N * (PFA.^(-1/N)-1);

caNME = sum(magSqGaussNoise)/N;

pd = sum(signalData > T * caNME)/numTrials;
