function [clutterEdgeData] = GenerateClutterEdge(noisedB, N, numTrials, numInterferers, CNR)
    if(numInterferers==0)
        clutterEdgeData = GenerateComplexGaussianNoise(noisedB, N, numTrials);
    else
        clutterEdgeData = GenerateComplexGaussianNoise(noisedB, N, numTrials, numInterferers, CNR);
    end
end
