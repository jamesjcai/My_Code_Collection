function [coords, loss_history] = dmercator(adjMatrix, D, maxIter, lr)
% D-Mercator network embedding in (D+1)-dimensional hyperbolic space
%
% INPUTS:
%   adjMatrix : NxN adjacency matrix (binary or weighted)
%   D         : similarity dimension (integer >= 1)
%   maxIter   : maximum number of optimization iterations
%   lr        : learning rate for gradient descent
%
% OUTPUTS:
%   coords        : Nx(D+1) matrix of node coordinates in hyperbolic space
%   loss_history  : vector of loss values over iterations
%
% Reference: D-Mercator method for multidimensional hyperbolic embedding

    % Validate inputs
    if nargin < 4
        error('Usage: dmercator(adjMatrix, D, maxIter, lr)');
    end
    if size(adjMatrix,1) ~= size(adjMatrix,2)
        error('Adjacency matrix must be square.');
    end
    if D < 1 || floor(D) ~= D
        error('D must be a positive integer.');
    end

    N = size(adjMatrix, 1);

    % Initialize coordinates: radial + angular (D-sphere)
    r = rand(N, 1) * 2; % radial coordinates
    theta = rand(N, D) * 2 * pi; % angular coordinates

    % Convert to Cartesian coordinates in (D+1)-dim hyperbolic space
    coords = hyperspherical_to_cartesian(r, theta);

    loss_history = zeros(maxIter, 1);

    % Optimization loop
    for iter = 1:maxIter
        % Compute pairwise hyperbolic distances
        distMat = hyperbolic_distance(coords);

        % Compute loss (negative log-likelihood)
        loss = embedding_loss(adjMatrix, distMat);
        loss_history(iter) = loss;

        % Compute gradients (simplified placeholder)
        grad = compute_gradients(adjMatrix, coords, distMat);

        % Update coordinates
        coords = coords - lr * grad;

        % Optional: reproject to valid hyperbolic space
        coords = reproject_to_hyperbolic(coords);

        % Display progress
        if mod(iter, 50) == 0
            fprintf('Iter %d, Loss: %.6f\n', iter, loss);
        end
    end
end

%% --- Helper Functions ---

function coords = hyperspherical_to_cartesian(r, theta)
    % Convert radial + angular coordinates to Cartesian in D+1 space
    N = length(r);
    D = size(theta, 2);
    coords = zeros(N, D+1);
    for i = 1:N
        coords(i,1) = cosh(r(i)); % time-like coordinate
        prod_sin = 1;
        for j = 1:D
            coords(i,j+1) = sinh(r(i)) * prod_sin * cos(theta(i,j));
            prod_sin = prod_sin * sin(theta(i,j));
        end
    end
end

function distMat = hyperbolic_distance(coords)
    % Compute pairwise hyperbolic distances using Minkowski metric
    N = size(coords,1);
    distMat = zeros(N);
    for i = 1:N
        for j = i+1:N
            minkowski_dot = coords(i,1)*coords(j,1) - sum(coords(i,2:end).*coords(j,2:end));
            dist = acosh(max(1, minkowski_dot));
            distMat(i,j) = dist;
            distMat(j,i) = dist;
        end
    end
end

function loss = embedding_loss(adjMatrix, distMat)
    % Simplified loss: logistic likelihood
    gamma = 1; % decay parameter
    P = 1 ./ (1 + exp((distMat - 1)/gamma));
    loss = -sum(sum(adjMatrix .* log(P + eps) + (1 - adjMatrix) .* log(1 - P + eps)));
end

function grad = compute_gradients(adjMatrix, coords, distMat)
    % Placeholder gradient computation (needs full derivation for real use)
    grad = randn(size(coords)) * 0.01; % random small perturbation
end

function coords = reproject_to_hyperbolic(coords)
    % Ensure coords remain on hyperboloid model
    for i = 1:size(coords,1)
        space_norm = norm(coords(i,2:end));
        coords(i,1) = sqrt(1 + space_norm^2);
    end
end
