function [estPFA] = ValidatePFA(noisedB, numTrials, NME, T)
% noisedB - noise power in decibels
% N - reference window size
% numTrials - number of trials
% NME - noise mean estimate
% T - threshold multiplier
% complexNoise - generated complex Gaussian noise data
% estPFA - estimated PFA

N = 1;
noisePower = 10^(noisedB/10);
complexNoise = sqrt(noisePower/2) * (randn(N,numTrials) + sqrt(-1) * randn(N,numTrials));
magSqGaussNoise = abs(complexNoise).^2;
estPFA = sum(magSqGaussNoise > NME * T)/numTrials;
end
