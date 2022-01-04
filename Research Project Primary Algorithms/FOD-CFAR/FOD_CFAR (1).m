function [pd, numCellsCensored] = FOD_CFAR(numTrials, PFA, N, SNR_dB, numInterferers)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% p - p lowest samples in the reference window (Number of homogeneous
% samples)
% a - threshold multiplier
% Ba - threshold
% SNR_dB - vector of SNR values in dB
% numInterferers - number of interferers in the reference window

% pd - probability of detection


numCellsCensored = zeros(numTrials);
%See FOD Paper for alpha value calculation (P error in paper in 0.2 and they use N=36 and p=24)
p = 24;
a = 5.9;
thresholdRange = 0.01:0.01:10;
noisedB = 0;

%Generating Signal Data that will be used to evaluate performance
numSNR = length(SNR_dB);
signalData = GenerateSignal(numTrials, SNR_dB);
pd = zeros(numSNR,1);
TFOD = zeros(1,N);
%Calculating Threshold for every number of homogeneous samples
%in the reference window
for k = 1:N
    TFOD(k) = FOD_Threshold(thresholdRange, k, N, PFA);
end
%Actual Algorithm
for n = 1:numSNR
    %Checking if interferer argument is present; If it is, it can be taken
    %into account when generating the noise data
    if (nargin == 4)
        magSqGaussNoise = sort(GenerateComplexGaussianNoise(noisedB, N, numTrials));
    elseif (nargin == 5)
        magSqGaussNoise = sort(GenerateComplexGaussianNoise(noisedB, N, numTrials, numInterferers, SNR_dB(n)));
    end
    %Calculating difference threshold and storing all differences between
    %consecutive samples in each trial
    Ba = a * std(magSqGaussNoise(1:p,:));
    differences = diff(magSqGaussNoise);

    %Splitting Homongeneous and Non-Homogeneous and Getting CA-CFAR NME for the
    %Homogeneous portion of each trial - Get an NME for every trial
    caNME = zeros(1,numTrials);
    Thresh = zeros(1, numTrials);
    for t = 1:numTrials
        numHomogeneousSamples = find(differences(:,t)>Ba(t),1);
        numCellsCensored(t) = N - numHomogeneousSamples;
        if(~isempty(numHomogeneousSamples))
            caNME(t) = sum(magSqGaussNoise(1:numHomogeneousSamples,t));
        %If numHomogeneousSamples has no value, the entire window is homogeneous, as there was no point in the window where the 
        %different between consecutive samples exceeded the different threshold (Ba), so you
        %can use the entire window for the NME calculation
        else
            numHomogeneousSamples = N;
            caNME(t) = sum(magSqGaussNoise(1:numHomogeneousSamples,t));
        end
        Thresh(t) = TFOD(numHomogeneousSamples);
    end

    %Threshold Multiplier Selection - Select the threshold that was
    %calculated using the number of homogeneous samples for that trial
    
    %Evaluating Performance
    pd(n) = sum(signalData(n,:) > (Thresh .* caNME))/numTrials;
end