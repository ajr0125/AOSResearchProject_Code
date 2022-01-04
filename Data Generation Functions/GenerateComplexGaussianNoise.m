function [magSqGaussNoise] = GenerateComplexGaussianNoise(noisedB, N, numTrials)
% noisedB - noise power in decibels
% N - reference window size
% numTrials - number of trials

% complexNoise - generated complex Gaussian noise data
noisePower = 10^(noisedB/10);
complexNoise = sqrt(noisePower/2) * (randn(N,numTrials) + sqrt(-1) * randn(N,numTrials));
magSqGaussNoise = abs(complexNoise).^2;

end
