function [T] = FOD_Threshold(possT, k, N, PFA)
% possT - possible Threshold values based on the calculated PFA value
% k - number of homogeneous samples in the reference window
% N - reference window size
% PFA - probability of false alarm goal

% T - calculated threshold multiplier
    possPFA = zeros(length(possT),1);
    %Condition for formula - k has to be greater than N/2, or more than
    %half the reference window
    if (k<= N/2)
        k = (N/2) + 1;
    end
    for m = 1:length(possT)
        j = 1:k;

        %FOD-CFAR Threshold Formula
        term = (possT(m) + (N-j+1)./(k-j+1)).^-1;
        possPFA(m) = nchoosek(N,k) * prod(term);
    end
    T = interp1(possPFA, possT, PFA);
end
