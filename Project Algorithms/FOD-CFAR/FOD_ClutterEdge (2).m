function [clutterEdgePFA, numCellsCensored] = FOD_ClutterEdge(numTrials, PFA, N, CNR)
% N - reference window size
% numTrials - number of trials
% PFA - probability of false alarm
% p - p lowest samples in the reference window (Number of homogeneous
% samples)
% a - threshold multiplier
% Ba - threshold
% numInterferers - number of interferers in the reference window
% CNR - clutter to noise ratio

% pd - probability of detection


numCellsCensored = zeros(numTrials);
%See FOD Paper for alpha value calculation (P error in paper in 0.2 and they use N=36 and p=24)
p = 24;
a = 5.9;
thresholdRange = 0.01:0.01:10;
noisedB = 0;
noisePower = 10^(noisedB/10);
clutterPower = 10^(CNR/10);

%Generating Data
clutterData = sort(exprnd(clutterPower, N, numTrials));
noiseData = sort(exprnd(noisePower, N, numTrials));
clutterNoiseData = noiseData;
iCUTClutter = GenerateSignal_ClutterEdge(numTrials,CNR);
iCUTSignal = GenerateSignal_ClutterEdge(numTrials, noisedB);
%Matrices for PFA and Thresholds at each number of interferers
clutterEdgePFA = zeros(N+1,1);
TFOD = zeros(1,N);
%Calculating Threshold for every number of homogeneous samples
%in the reference window
for k = 1:N
    TFOD(k) = FOD_Threshold(thresholdRange, k, N, PFA);
end
%Actual Algorithm
for numInterferers = 0:N
    %Generating Signal (CUT) Accordingly - If the majority of the window is
    %clutter, the signal will also be clutter, else it will be normal noise
    %(signal)
    if(numInterferers>N/2)
        iCUT = iCUTClutter;
    else
        iCUT = iCUTSignal;
    end
    %Assigning Interferers to Data
    clutterNoiseData(1:numInterferers,:) = clutterData(1:numInterferers,:);

    %Calculating difference threshold and storing all differences between
    %consecutive samples in each trial
    Ba = a * std(clutterNoiseData(1:p,:));
    differences = diff(clutterNoiseData);

    %Splitting Homongeneous and Non-Homogeneous and Getting CA-CFAR NME for the
    %Homogeneous portion of each trial - Get an NME for every trial
    caNME = zeros(1,numTrials);
    Thresh = zeros(1, numTrials);
    for t = 1:numTrials
        numHomogeneousSamples = find(differences(:,t)>Ba(t),1);
        numCellsCensored(t) = N - numHomogeneousSamples;
        if(~isempty(numHomogeneousSamples))
            caNME(t) = sum(clutterNoiseData(1:numHomogeneousSamples,t));
        %If numHomogeneousSamples has no value, the entire window is homogeneous, as there was no point in the window where the 
        %different between consecutive samples exceeded the different threshold (Ba), so you
        %can use the entire window for the NME calculation
        else
            numHomogeneousSamples = N;
            caNME(t) = sum(clutterNoiseData(1:numHomogeneousSamples,t));
        end
        %Threshold Multiplier Selection - Select the threshold that was
        %calculated using the number of homogeneous samples for that trial
        Thresh(t) = TFOD(numHomogeneousSamples);
    end
    %Evaluating PFA
    clutterEdgePFA(numInterferers+1) = sum(iCUT > (Thresh .* caNME))/numTrials;
end