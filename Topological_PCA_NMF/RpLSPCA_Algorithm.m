Here is the MATLAB translation of the provided Python code:

```matlab
function gaussian_kernel = gaussian_kernel(dist, t)
%gaussian kernel function for weighted edges
gaussian_kernel = exp(-(dist.^2 / t));
end

function dist_mat = Eu_dis(x)
%Calculate the distance among each row of x
%x: N X D, N: number of samples, D: Dimension of the feature
%Return: N X N distance matrix
x = double(x);
aa = sum(x.*x, 2);
ab = x * x';
dist_mat = aa + aa' - 2 * ab;
dist_mat(dist_mat < 0) = 0;
dist_mat = sqrt(dist_mat);
dist_mat = max(dist_mat, dist_mat');
dist_mat = double(dist_mat);
end

function W_L = cal_weighted_adj(data, n_neighbors, t)
%Calculate weighted adjacency matrix based on KNN
%For each row of X, put an edge between nodes i and j
%If nodes are among the n_neighbors nearest neighbors of each other
%according to Euclidean distance
dist = Eu_dis(data);
n = size(dist, 1);
gk_dist = gaussian_kernel(dist, t);
W_L = zeros(n, n);
for i = 1:n
    index_L = sort(dist(i, :));
    index_L = index_L(2:n_neighbors+1);
    len_index_L = length(index_L);
    for j = 1:len_index_L
        W_L(i, index_L(j)) = gk_dist(i, index_L(j)); %weighted edges
    end
end
W_L = max(W_L, W_L');
end

function W_L = cal_unweighted_adj(data, n_neighbors)
%Calculate unweighted adjacency matrix based on KNN
dist = Eu_dis(data);
n = size(dist, 1);
W_L = zeros(n, n);
for i = 1:n
    index_L = sort(dist(i, :));
    index_L = index_L(2:n_neighbors+1);
    len_index_L = length(index_L);
    for j = 1:len_index_L
        W_L(i, index_L(j)) = 1; %edges not weighted
    end
end
W_L = max(W_L, W_L');
end

function L = cal_laplace(adj)
N = size(adj, 1);
D = zeros(N, N);
for i = 1:N
    D(i, i) = sum(adj(i, :)); %Degree Matrix
end
L = D - adj; %Laplacian
end

function U = RpLSPCA_Algorithm(xMat, laplace, beta, gamma, k, n)
%Optimization Algorithm of RgLSPCA / RpLSPCA
%Solve approximately via ADMM
%Need to compute optimal principal directions matrix U
%Projected Data matrix Q
%Error term matrix E = X - UQ^T
%Z matrix used to solve Q (see supplementary information)
%Inputs are data matrix X, laplacian term, scale parameters,
%number of reduced dimensions, number of original dimensions

%Initialize thresholds, matrices
obj1 = 0;
obj2 = 0;
thresh = 1e-50;
V = eye(n);
vMat = double(V); %Auxiliary matrix to optimize L2,1 norm
E = ones(size(xMat));
E = double(E); %Error term X - UQ^T
C = ones(size(xMat));
C = double(C); %Lagrangian Multiplier
laplace = double(laplace); %Laplacian
miu = 1; %Penalty Term
for m = 1:30
    Z = (-(miu/2) * ((E - xMat + C/miu)' * (E - xMat + C/miu))) + beta * vMat + gamma * laplace;
    %cal Q (Projected Data Matrix)
    [Z_eigVects, Z_eigVals] = eig(double(Z));
    [~, eigValIndice] = sort(diag(Z_eigVals));
    n_eigValIndice = eigValIndice(1:k);
    n_Z_eigVect = Z_eigVects(:, n_eigValIndice);
    %Optimal Q given by eigenvectors corresponding
    %to smallest k eigenvectors
    Q = n_Z_eigVect;
    %cal V
    q = vecnorm(Q, 2, 2);
    qq = 1 ./ (q * 2);
    VV = diag(qq);
    vMat = double(VV);
    qMat = double(Q);
    %cal U (Principal Directions)
    U = (xMat - E - C/miu) * qMat;
    %cal P (intermediate step)
    P = xMat - U * qMat' - C/miu;
    %cal E (Error Term)
    for i = 1:size(E, 2)
        E(:, i) = (max((1 - 1 / (miu * vecnorm(P(:, i)))), 0)) * P(:, i);
    end
    %update C
    C = C + miu * (E - xMat + U * qMat');
    %update miu
    miu = 1.2 * miu;

    obj1 = vecnorm(qMat);
    if m > 1
        diff = obj2 - obj1;
        if diff < thresh
            break; %end iterations if error within accepted threshold
        end
    end
    obj2 = obj1;
end
end

function P_L = cal_persistent_laplace(W_L, zetas)
n = size(W_L, 1);
W_L(logical(eye(n))) = 0;

L = cal_laplace(W_L);
L(logical(eye(n))) = 1e8; %Make sure diagonal is excluded from maximal and minimal value consideration
min_l = min(L(L ~= 0)); %Establish Min Value
L(logical(eye(n))) = -1e8;
max_l = max(L(L ~= 0)); %Establish Max Value

d = max_l - min_l;

L = cal_laplace(W_L);
PL = zeros(8, n, n);
for k = 1:7
    PL(k, :, :) = double(L < (k/7*d + min_l));
    PL(k, logical(eye(n))) = 0;
    PL(k, :, :) = cal_laplace(PL(k, :, :));
end

P_L = sum(zetas(:, :, :) .* PL, 1);
end

function Y = RpLSPCA_cal_projections(X_data, beta1, gamma1, k_d, zetas)
n = size(X_data, 1);
dist = Eu_dis(X_data);
max_dist = max(max(dist));
W_L = cal_weighted_adj(X_data, 15, max_dist^(2));
A = W_L;
PL = cal_persistent_laplace(A, zetas);
Y = RpLSPCA_Algorithm(X_data', PL, beta1, gamma1, k_d, n);
end

function Persistent_Laplacian = cal_persistence_KNN(data, n_filtrations, zetas)
n = size(data, 1);
%Consider n neighbors and reduce by 2 neighbors at
%each iteration of filtration down to 1 nearest neighbor
%(p filtrations)
num_neighbors_list = (1:2:n_filtrations)';
num_filtrations = length(num_neighbors_list);

PL = zeros(num_filtrations, n, n);
zetas = double(zetas);

for idx = 1:num_filtrations
    num_neighbors = num_neighbors_list(idx);
    A = cal_unwe