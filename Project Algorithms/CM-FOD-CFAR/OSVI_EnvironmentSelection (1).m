function [noiseData, TSelect] = OSVI_EnvironmentSelection(N, magSqGaussNoise, KVI, KMR)
% N - reference window size
% PFA - probability of false alarm
% magSqGaussNoise - noise data
% signalData - signal data
% k - ordered index value

% pd - probability of detection


%Calculating Variables for VI and MR
leadWindow = magSqGaussNoise(1:(N/2));
lagWindow = magSqGaussNoise((N/2)+1:N);
lagSum = sum(lagWindow);
leadSum = sum(leadWindow);
lagVI = 1 + var(lagWindow)/(mean(lagWindow)^2);
leadVI = 1 + var(leadWindow)/(mean(leadWindow)^2);
MR = leadSum/lagSum;

%Comparing VI and MR to their thresholds
leadVariable = leadVI>KVI;
lagVariable = lagVI>KVI;
diffMean = (KMR^-1>MR) | (MR>KMR);
oneVariable = xor(leadVariable, lagVariable);
bothVariable = leadVariable & lagVariable;
noneVariable = (~leadVariable) & (~lagVariable);

%Homogeneous Case - Same Mean, no windows variable (OS-CFAR on both
%windows) - Use TFull (0) and the entire window
 if ~diffMean && noneVariable
    noiseData = magSqGaussNoise;
    TSelect = 0;
 end

%Clutter Case - Different mean, no windows variable (OSGO-CFAR) - Use
%TClutter (1) and the window with a greater power sum
if diffMean && noneVariable
    if leadSum>lagSum
        noiseData = sort(leadWindow);
    else
        noiseData = sort(lagWindow);
    end
    TSelect = 1;
end


%OS with one Window - Leading window variable, lagging window not variable
%(OS-CFAR on Lagging Window) - Use THalf (2) and the lagging window
if oneVariable && leadVariable
    noiseData = lagWindow;
    TSelect = 2;
end


%OS with one Window - Leading window not variable, lagging window variable
%(OS-CFAR on Leading Window) - Use THalf (2) and the leading window
if oneVariable && lagVariable
   noiseData = leadWindow;
   TSelect = 2;
end

%%Both windows variable - OSSO - Use THalf (2) and the window with the
%%smaller power sum
if bothVariable
    if leadSum<lagSum
        noiseData = leadWindow;
    else
        noiseData = lagWindow;
    end
    TSelect = 2;
end




