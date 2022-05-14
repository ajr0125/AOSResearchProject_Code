function [pd] = SO_CFAR(numTrials, PFA, N, magSqGaussNoise, signalData)
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
    possPFA(n) = 2 * (2 + (possT(n)/(N/2))).^(-N/2) * eqSum;
end
T = interp1(possPFA, possT, PFA);
TSO = T/(N/2);

lagSum = sum(magSqGaussNoise(1:(N/2),:));
leadSum = sum(magSqGaussNoise((N/2)+1:N,:));
soNME = min(lagSum,leadSum);

pd = sum(signalData > TSO * soNME)/numTrials;