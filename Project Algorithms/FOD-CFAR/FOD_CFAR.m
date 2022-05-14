function [pd, numCellsCensored] = FOD_CFAR(numTrials, PFA, N, magSqGaussNoise, signalData)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% p - p lowest samples in the reference window (Number of homogeneous
% samples)
% a - threshold multiplier
% Ba - threshold
% magSqGaussNoise - noise data
% signalData - signal data

% pd - probability of detection


numCellsCensored = zeros(1,numTrials);
%See FOD Paper for alpha value calculation (P error in paper in 0.2 and they use N=36 and p=24)
p = 24;
a = 5.9;
thresholdRange = 0.01:0.01:10;

TFOD = zeros(1,N);
%Calculating Threshold for every number of homogeneous samples
%in the reference window
for k = 1:N
    TFOD(k) = FOD_Threshold(thresholdRange, k, N, PFA);
end
%Actual Algorithm
%Sorting Noise Data
magSqGaussNoise = sort(magSqGaussNoise);
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
    if(~isempty(numHomogeneousSamples))
        numCellsCensored(t) = N - numHomogeneousSamples;
        caNME(t) = sum(magSqGaussNoise(1:numHomogeneousSamples,t));
    %If numHomogeneousSamples has no value, the entire window is homogeneous, as there was no point in the window where the 
    %different between consecutive samples exceeded the different threshold (Ba), so you
    %can use the entire window for the NME calculation
    else
        numHomogeneousSamples = N;
        numCellsCensored(t) = N - numHomogeneousSamples;
        caNME(t) = sum(magSqGaussNoise(1:numHomogeneousSamples,t));
    end
    %Threshold Multiplier Selection - Select the threshold that was
    %calculated using the number of homogeneous samples for that trial
    Thresh(t) = TFOD(numHomogeneousSamples);
end

%Evaluating Performance
pd = sum(signalData > (Thresh .* caNME))/numTrials;
