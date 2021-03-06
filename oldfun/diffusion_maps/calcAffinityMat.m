function [K, nnData] = calcAffinityMat(X, configParams)
% calculate affinity matrix of data matrix X, size N by M 
% N - length of feature vector, M - number of vectors
% affinity matrix is calculated for kNN nearest neighbors, resulting in
% sparse matrix. This saves on runtime.
% the scale for the affinity matrix can be set using auto-tuning
%
% Gal Mishne

dParams.dist_type = 'euclidean';
dParams.kNN = 4;
dParams.self_tune = 0;
dParams.verbose = false;
configParams = setParams(dParams, configParams);

[~,M] = size(X);
configParams.kNN = min(configParams.kNN,M);

%% affinity matrix
tic
[Dis, Inds] = pdist2(X',X',configParams.dist_type,'Smallest',configParams.kNN);
toc
            
nnData.nnDist = Dis;
nnData.nnInds = Inds;

% Total number of entries in the W matrix
numNonZeros = configParams.kNN * M;
% Using sparse matrix for affinity matrix
ind = 1;

% calc the sparse row and column indices
rowInds = repmat((1:M),configParams.kNN,1);
rowInds = rowInds(:);
colInds = double(nnData.nnInds(:));
vals    = nnData.nnDist(:);

autotuneVals = zeros(numNonZeros,1);
if configParams.self_tune
    nnAutotune = min(configParams.self_tune,size(nnData.nnDist,2));
    sigmaKvec = (nnData.nnDist(nnAutotune,:));
end
if configParams.self_tune
    for i = 1:M
        autotuneVals(ind : ind + configParams.kNN - 1) = sigmaKvec(i) * sigmaKvec(nnData.nnInds(:,i));
        ind = ind + configParams.kNN;
    end
end
% sparse affinity matrix of dataset X
nnData.dstsX = sparse(rowInds, colInds, vals, M, M);
nnData.indsXX = sub2ind([M, M], rowInds, colInds);


% setting local scale for each sample, 
% following "Self-tuning Clustering", Zelink-Manor and Perona
if configParams.self_tune
    vals = exp(-vals.^2./(autotuneVals+eps));
    K = sparse(rowInds, colInds, vals, M, M);
else
    sig = median(vals);
    vals = exp(-vals.^2/sig^2);
    K = sparse(rowInds, colInds, vals, M, M);
end
K = (K + K')/2;
return
