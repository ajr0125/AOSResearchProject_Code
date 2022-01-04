function [pd] = GO_CFAR(numTrials, PFA, N, SNR_dB)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% T - threshold multiplier
% SNR_dB - vector of SNR values in dB

% pd - probability of detection
possT = 0.1:0.1:50;
possPFA = zeros(length(possT),1);
for n = 1:length(possT)
    eqSum = 0;
    for k = 0:(N/2)-1
        eqSum = eqSum + nchoosek((N/2)-1+k, k) * (2 + (possT(n)/(N/2))).^(-k);
    end
    possPFA(n) = 2 * ((1 + (possT(n)/(N/2))).^(-N/2) - (2 + (possT(n)/(N/2))).^(-N/2) * eqSum);
end
T = interp1(possPFA, possT, PFA);
TGO = T/(N/2);
noisedB = 0;
magSqGaussNoise = GenerateComplexGaussianNoise(noisedB, N, numTrials);

lagSum = sum(magSqGaussNoise(1:(N/2),:));
leadSum = sum(magSqGaussNoise((N/2)+1:N,:));
goNME = min(lagSum,leadSum);

numSNR = length(SNR_dB);
signalData = GenerateSignal(numTrials, SNR_dB);
pd = zeros(numSNR,1);
for n = 1:numSNR
    pd(n) = sum(signalData(n,:) > TGO * goNME)/numTrials;
end
