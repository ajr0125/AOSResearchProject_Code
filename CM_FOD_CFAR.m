%Clutter Mitigated First Order Difference CFAR (CM-FOD-CFAR)
function [pd] = CM_FOD_CFAR(PFA, magSqGaussNoise, signalData, KVI, KMR, kHalf, p, a)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% magSqGaussNoise - noise data
% signalData - signal data
% KVI - OSVI-CFAR hypothesis testing threshold for variability index
% KMR - OSVI-CFAR hypothesis testing threshold for mean ratio
% kHalf - OSVI-CFAR: kth ordered index value for the clutter case
% p - FOD-CFAR: p lowest samples in the reference window (Number of homogeneous
% samples)
% a - FOD-CFAR: threshold multiplier

% pd - probability of detection

%Establishing Data Dimensions and Array for Probability of Detection
[N,numTrials] = size(magSqGaussNoise);
pdArray = zeros(1, numTrials);

%Threshold Calculation for FOD-CFAR: Array of Thresholds for every possible
%number of interferers
thresholdRange = 0.01:0.01:10;
TFOD = zeros(1,N);
for k = 1:N
    TFOD(k) = FOD_Threshold(thresholdRange, k, N, PFA);
end

%Threshold Calculation for OSVI-CFAR
%Clutter Case - take the greatest value in one window as NME (OSGO)
%The clutter case is the case the FOD-CFAR is weakest in, thus we use the
%OSVI-CFAR instead
TClutter = CalculateOSThreshold(PFA, N/2, N/2);
clutterCount = 0;
tic
for t = 1:numTrials
    %Identify the noise data environment
    [noiseData, TSelect] = OSVI_EnvironmentSelection(N, magSqGaussNoise(:,t), KVI, KMR);

    %If in the clutter case, use the OSVI-CFAR clutter threshold and kHalf
    %to calculate the noise mean estimate and the threshold
    if TSelect==1
        pdArray(t) = signalData(:,t)>(TClutter * noiseData(kHalf));
        clutterCount = clutterCount+1;

    %If not in the clutter case, use regular FOD-CFAR
    else
        pdArray(t) = FOD_OSDVI_Version(p, a, TFOD, magSqGaussNoise(:,t), signalData(:,t));
    end
end
%Evaluate Performance
pd = sum(pdArray)/numTrials;
disp(['Clutter Count: ' num2str(clutterCount)])
toc



