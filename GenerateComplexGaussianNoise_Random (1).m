function [magSqGaussNoise] = GenerateComplexGaussianNoise_Random(noisedB, N, numTrials, numInterferers, INR, isSplit)
% noisedB - noise power in decibels
% N - reference window size
% numTrials - number of trials
% numInterferers - number of interferers
% INR - mean power of interferers, if present
% isSplit - indicates whether the interferers should be in both windows or
% just one

% complexNoise - generated complex Gaussian noise data

% magSqGaussNoise - noise data
if (nargin == 3)
    noisePower = 10^(noisedB/10);
    complexNoise = sqrt(noisePower/2) * (randn(N,numTrials) + sqrt(-1) * randn(N,numTrials));
    magSqGaussNoise = abs(complexNoise).^2;
elseif (nargin == 5)
    interfererPower = 10.^(INR/10);

    %Generating Swerling 1 Model for Interferers (Varying around
    %interfererPower but the mean is equal to interfererPower)
    intData = exprnd(interfererPower,numInterferers,numTrials);

    noisePower = 10^(noisedB/10);
    complexNoise = sqrt(noisePower/2) * (randn(N,numTrials) + sqrt(-1) * randn(N,numTrials));
    indices = ceil(rand(1,numInterferers) * N/2);
    while(unique(indices)~=numInterferers)
        indices = ceil(rand(1,numInterferers) * N/2);
    end
    complexNoise(indices,:) = sqrt(intData/2) + sqrt(-1) * sqrt(intData/2);
    magSqGaussNoise = abs(complexNoise).^2;
    %magSqGaussNoise(indices,:) = exprnd(interfererPower, numInterferers, numTrials);
elseif (nargin==6)
    if(rem(numInterferers,2)==1)
        perror('Number of Interferers must be even in split window case');
    end
    interfererPower = 10.^(INR/10);

    %Generating Swerling 1 Model for Interferers (Varying around
    %interfererPower but the mean is equal to interfererPower)
    intData = exprnd(interfererPower,numInterferers,numTrials);

    noisePower = 10^(noisedB/10);
    complexNoise = sqrt(noisePower/2) * (randn(N,numTrials) + sqrt(-1) * randn(N,numTrials));
    %magSqGaussNoise = exprnd(1,N,numTrials);
    indices = ceil(rand(1,numInterferers/2) * N/2);
    while(unique(indices)~=numInterferers/2)
        indices = ceil(rand(1,numInterferers/2) * N/2);
    end
    complexNoise(indices,:) = sqrt(intData(1:numInterferers/2,:)/2) + sqrt(-1) * sqrt(intData(1:numInterferers/2,:)/2);
    complexNoise(indices+N/2,:) = sqrt(intData((numInterferers/2)+1:numInterferers,:)/2) + sqrt(-1) * sqrt(intData((numInterferers/2)+1:numInterferers,:)/2);
    magSqGaussNoise = abs(complexNoise).^2;
    %magSqGaussNoise(indices,:) = intData(1:numInterferers/2,:);
    %magSqGaussNoise(indices+N/2,:) = intData((numInterferers/2)+1:numInterferers,:);

else
    disp('Error: Wrong number of arguments');
end
end