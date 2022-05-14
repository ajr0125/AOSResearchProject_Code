function [pd] = OSSO_CFAR(numTrials, PFA, N, magSqGaussNoise, signalData, k)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% T - threshold multiplier
% SNR_dB - vector of SNR values in dB
% k - ordered index value

% pd - probability of detection
possT = 0.1:0.1:50;
possPFA = zeros(length(possT),1);
for n = 1:length(possT)
    possPFA(n) = prod((1 + (possT(n))./(N + 1 - (1:k))).^-1);
    %possPFA(n) = k * nchoosek(N,k) * factorial(k-1) * factorial(floor(possT(n) + N - k)-1) * (possT(n) + N - k)/(factorial(floor(possT(n) + N)-1) * (possT(n) + N));
end
T = interp1(possPFA, possT, PFA);
magSqGaussNoise = sort(magSqGaussNoise);

%Storing all data in window variables and summing up the windows
lagWindow = magSqGaussNoise(1:(N/2),:);
leadWindow = magSqGaussNoise((N/2)+1:N,:);
lagSum = sum(magSqGaussNoise(1:(N/2),:));
leadSum = sum(magSqGaussNoise((N/2)+1:N,:));

%Taking all lead window k values, then replacing the ones in the trials
%where the lagging window sum is less than the leading window sum
ossoNME = leadWindow(k,:);
idx = lagSum<leadSum;
ossoNME(idx) = lagWindow(k,idx);

pd = sum(signalData > T * ossoNME)/numTrials;