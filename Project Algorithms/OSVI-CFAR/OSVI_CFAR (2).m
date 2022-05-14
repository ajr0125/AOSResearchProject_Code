function [pd] = OSVI_CFAR(PFA, magSqGaussNoise, signalData, KVI, KMR, kFull, kHalf)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% magSqGaussNoise - noise data
% signalData - signal data
% k - ordered index value

% pd - probability of detection
[N,numTrials] = size(magSqGaussNoise);
%Threshold Calculation
%Case for Taking only One Window
THalf = CalculateOSThreshold(PFA, N/2, kHalf);
%Case for taking both windows
TFull = CalculateOSThreshold(PFA, N, kFull);
%Case for taking greatest value in one window as NME (OSGO)
TClutter = CalculateOSThreshold(PFA, N/2, N/2);


%Calculating Variables for VI and MR
leadWindow = magSqGaussNoise(1:(N/2),:);
lagWindow = magSqGaussNoise((N/2)+1:N,:);
lagSum = sum(lagWindow);
leadSum = sum(leadWindow);
lagVI = 1 + var(lagWindow)./(mean(lagWindow).^2);
leadVI = 1 + var(leadWindow)./(mean(leadWindow).^2);
MR = leadSum./lagSum;

%Comparing VI and MR to their thresholds
leadVariable = leadVI>KVI;
lagVariable = lagVI>KVI;
diffMean = (KMR^-1>MR) | (MR>KMR);
oneVariable = xor(leadVariable, lagVariable);
bothVariable = leadVariable & lagVariable;
noneVariable = (~leadVariable) & (~lagVariable);

%Sorted Windows
lagWindowSorted = sort(lagWindow);
leadWindowSorted = sort(leadWindow);
fullWindowSorted = sort(magSqGaussNoise);

%All Noise Mean Estimates
nmeHomogeneous = fullWindowSorted(kFull,:);
lagNME = lagWindowSorted(kHalf,:);
leadNME = leadWindowSorted(kHalf,:);
lagNMEGO = lagWindowSorted(N/2,:);
leadNMEGO = leadWindowSorted(N/2,:);

%Homogeneous Case - Same Mean, no windows variable (OS-CFAR on both
%windows)
idxHomogeneous = ~diffMean & noneVariable;
Threshold(idxHomogeneous) = TFull * nmeHomogeneous(idxHomogeneous); 

%Clutter Case - Different mean, no windows variable (OSGO-CFAR)
idxClutter = diffMean & noneVariable;
Threshold(idxClutter) = (((lagSum(idxClutter)>leadSum(idxClutter)) .* lagNMEGO(idxClutter)) ...
    + ((lagSum(idxClutter)<leadSum(idxClutter)) .* leadNMEGO(idxClutter))) * TClutter;


%OS with one Window - Leading window variable, lagging window not variable
%(OS-CFAR on Lagging Window)
idxLag = oneVariable & leadVariable;
Threshold(idxLag) = THalf * lagNME(idxLag);


%OS with one Window - Leading window not variable, lagging window variable
%(OS-CFAR on Leading Window)
idxLead = oneVariable & lagVariable;
Threshold(idxLead) = THalf * leadNME(idxLead);

%%Both windows variable - OSSO
idxVariable = bothVariable;
Threshold(idxVariable) = (((lagSum(idxVariable)>leadSum(idxVariable)) .* leadNME(idxVariable)) ...
    + ((lagSum(idxVariable)<leadSum(idxVariable)) .* lagNME(idxVariable))) * THalf;

pd = sum(signalData>Threshold)/numTrials;



