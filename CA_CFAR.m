function [pd] = CA_CFAR(numTrials, PFA, N, SNR_dB)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% T - threshold multiplier
% SNR_dB - vector of SNR values in dB

% pd - probability of detection
T = N * (PFA.^(-1/N)-1);
noisedB = 0;
magSqGaussNoise = GenerateComplexGaussianNoise(noisedB, N, numTrials);

caNME = sum(magSqGaussNoise)/N;

numSNR = length(SNR_dB);
signalData = GenerateSignal(numTrials, SNR_dB);
pd = zeros(numSNR,1);
for n = 1:numSNR
    pd(n) = sum(signalData(n,:) > T * caNME)/numTrials;
end