function [pd] = CMLD_CFAR(numTrials, PFA, N, SNR_dB)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% T - threshold multiplier
% SNR_dB - vector of SNR values in dB

% pd - probability of detection
j = input("Input a the number of interferers, ranging from 1 to " + N + ": ");
k = N - j;
noisedB = 0;
magSqGaussNoise = GenerateComplexGaussianNoise(noisedB, N, numTrials);
sortedmagSqGaussNoise = sort(magSqGaussNoise);
a = PFA.^(-1/k)-1;
summation = zeros(1, numTrials);
for i = 1:k-1
    summation = summation + sortedmagSqGaussNoise(i,:);
end
T = a * (sortedmagSqGaussNoise(12,:) * (N-k) + summation);

cmldNME = sum(sortedmagSqGaussNoise(1:k,:))/k;

numSNR = length(SNR_dB);
signalData = GenerateSignal(numTrials, SNR_dB);
pd = zeros(numSNR,1);
for n = 1:numSNR
    pd(n) = sum(signalData(n,:) > T .* cmldNME)/numTrials;
end
