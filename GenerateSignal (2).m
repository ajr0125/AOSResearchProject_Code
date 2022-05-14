function [powerData] = GenerateSignal(numTrials, SNR_dB)
% noisedB - noise power in decibels
% N - reference window size
% numTrials - number of trials
% NME - noise mean estimate
% T - threshold multiplier
% complexNoise - generated complex Gaussian noise data

% powerData - magnitude squared signal power data
noisePower = 1;
SNRCount = length(SNR_dB);
SNRPower = 10.^(SNR_dB/10);
powerData = zeros(SNRCount,numTrials);
noise = sqrt(noisePower/2) * (randn(1,numTrials) + sqrt(-1) * randn(1,numTrials));
for n = 1:SNRCount
    sw = exprnd(SNRPower(n), 1, numTrials);
    signal = sqrt(sw/2) + sqrt(-1) * sqrt(sw/2);
    powerData(n,:) = abs(signal+noise).^2;
end