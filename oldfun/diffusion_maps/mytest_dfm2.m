%% Demostration of Filter, Normalization and Batch Correction of Data in scGEApp
%% Read scRNA-seq data, X and Y
[X,genelistx]=sc_readfile('../../example_data/GSM3204304_P_P_Expr_999cells.csv');
[Y,genelisty]=sc_readfile('../../example_data/GSM3204305_P_N_Expr_999cells.csv');

%% Select genes with at least 3 cells having more than 5 reads per cell. 
[X,genelistx]=sc_selectg(X,genelistx,5,3);
[Y,genelisty]=sc_selectg(Y,genelisty,5,3);

%% Obtain gene intersection of X and Y
[genelist,i,j]=intersect(genelistx,genelisty,'stable');
X=X(i,:);
Y=Y(j,:);
% libsizex=sum(X);
% libsizey=sum(Y);
% X=X(:,libsizex>quantile(libsizex,0.3)&libsizex<quantile(libsizex,0.95));
% Y=Y(:,libsizey>quantile(libsizey,0.3)&libsizey<quantile(libsizey,0.95));
% clearvars -except X Y genelist

%%

D=pdist2(X',X','cosine');
dt=D-diag(diag(D));

N=999;
m=999;

% Changing these values will lead to different nonlinear embeddings
knn    = ceil(0.03*N); % each patch will only look at its knn nearest neighbors in R^d
sigma2 = 100; % determines strength of connection in graph... see below
%%
[srtdDt,srtdIdx] = sort(dt,'ascend');
dt               = srtdDt(1:knn+1,:);
nidx             = srtdIdx(1:knn+1,:);
cmap = jet(N);
%%

% compute weights
tempW  = exp(-dt.^2/sigma2); 
% build weight matrix
i = repmat(1:m,knn+1,1);
W = sparse(i(:),double(nidx(:)),tempW(:),m,m); 
W = max(W,W'); % for undirected graph.
% The original normalized graph Laplacian, non-corrected for density
ld = diag(sum(W,2).^(-1/2));
DO = ld*W*ld;
DO = max(DO,DO');%(DO + DO')/2;
% get eigenvectors
[v,d] = eigs(DO,10,'la');

figure;
scatter(v(:,1),v(:,2))
return;

eigVecIdx = nchoosek(2:4,2);
for i = 1:size(eigVecIdx,1)
    figure,scatter(v(:,eigVecIdx(i,1)),v(:,eigVecIdx(i,2)),20,cmap)
    title('Nonlinear embedding');
    xlabel(['\phi_',num2str(eigVecIdx(i,1))]);
    ylabel(['\phi_',num2str(eigVecIdx(i,2))]);
end





