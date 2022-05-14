function [performance] = AreaPerformance(SNR_dB, algorithmPD, PFA, signalPower)
% algorithmPD - Probability of detection for the tested algorithm
% PDOPT - probability of detection for the optimal curve

% area - Area under the curve of the tested algorithm
% totalArea - total area in the graph (length of x-axis * length of y-axis)
% areaPortion - portion of the graph area taken up by the curve
% optArea - area under the optimal curve
% optAreaPortion - portion of the graph area taken up by the optimal curve

% performance - ratio of algorithm's area portion to the optimal area
% portion


%Area under curve
area = trapz(SNR_dB, algorithmPD);

%Optimal area
Nopt = 1e6;
TOPT = (PFA^(-1/Nopt))-1;
PDOPT = (1 + (TOPT./(1+signalPower))).^(-Nopt);
optArea = trapz(SNR_dB, PDOPT);

performance = area/optArea;