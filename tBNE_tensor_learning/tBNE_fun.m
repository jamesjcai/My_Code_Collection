function [T, W] = tBNE_fun(X, Z, Y, k)
% t-BNE computes brain network embedding based on constrained tensor factorization
%
% INPUT
% X: brain networks stacked in a 3-way tensor
% Z: side information
% Y: label information
% k: rank of CP factorization
%
% OUTPUT
% T is the factor tensor containing
%   vertex factor matrix B = T{1} and
%   subject factor matrix S = T{3}
% W is the weight matrix
%
% Example: see tBNE_demo.m
%
% Reference:
% Bokai Cao, Lifang He, Xiaokai Wei, Mengqi Xing, Philip S. Yu, 
% Heide Klumpp and Alex D. Leow. t-BNE: Tensor-based Brain Network Embedding.
% In SDM 2017.
%
% Dependency:
% [1] Matlab tensor toolbox v 2.6
% Brett W. Bader, Tamara G. Kolda and others
% http://www.sandia.gov/~tgkolda/TensorToolbox 
% [2] A feasible method for optimization with orthogonality constraints
% Zaiwen Wen and Wotao Yin
% http://www.math.ucla.edu/~wotaoyin/papers/feasible_method_matrix_manifold.html

    %% set algorithm parameters
    printitn = 10;
    maxiter = 200;
    fitchangetol = 1e-4;

    alpha = 0.1; % weight for guidance
    beta = 0.1; % weight for classification loss
    gamma = 0.1; % weight for regularization

    u = 1e-6;
    umax = 1e6;
    rho = 1.15;

    opts.record = 0;
    opts.mxitr = 20;
    opts.xtol = 1e-5;
    opts.gtol = 1e-5;
    opts.ftol = 1e-8;

    %% compute statistics
    dim = size(X);
    normX = norm(X);
    numClass = size(Y, 2);
    m = dim(1);
    n = dim(3);
    l = size(Y, 1);
    D = [eye(l), zeros(l, n - l)];
    L = diag(sum(Z * Z')) - Z * Z';

    %% initialization
    B = randn(m, k);
    P = B;
    S = randn(n, k);
    S = orth(S);
    W = randn(k, numClass);
    U = zeros(m, k); % Lagrange multipliers

    %% main loop
    fit = 0;
    for iter = 1 : maxiter
        fitold = fit;
        % update B
        ete = (S' * S) .* (P' * P); % compute E'E
        b = 2 * ete + u * eye(k);
        c = 2 * mttkrp(X, {B, P, S}, 1) + u * P + U;
        B = c / b;

        % update P
        ftf = (S' * S) .* (B' * B); % compute F'F
        b = 2 * ftf + u * eye(k);
        c = 2 * mttkrp(X, {B, P, S}, 2) + u * B - U;
        P = c / b;

        % update U
        U = U + u * (P - B);

        % update u
        u = min(rho * u, umax);

        % update S
        tic;
        [S, out] = OptStiefelGBB(...
            S, @Sfun, opts, B, P, X, L, D, W, Y, alpha, beta);
        tsolve = toc;
        fprintf(...
            ['[S]: obj val %7.6e, cpu %f, #func eval %d, ', ...
            'itr %d, |ST*S-I| %3.2e\n'], ...
            out.fval, tsolve, out.nfe, out.itr, norm(S' * S - eye(k), 'fro'));

        % update W
        H = D * S;
        W = (H' * H + gamma * eye(k)) \ H' * Y;

        % compute the fit
        T = ktensor({B, P, S});
        normresidual = sqrt(normX ^ 2 + norm(T) ^ 2 - 2 * innerprod(X, T));
        fit = 1 - (normresidual / normX);
        fitchange = abs(fitold - fit);

        if mod(iter, printitn) == 0
            fprintf(' Iter %2d: fitdelta = %7.1e\n', iter, fitchange);
        end

        % check for convergence
        if (iter > 1) && (fitchange < fitchangetol)
            break;
        end
    end

    %% clean up final results
    T = arrange(T); % columns are normalized

    fprintf('factorization error %3.2e\n', fit);
end

function [F, G] = Sfun(S, B, P, X, L, D, W, Y, alpha, beta)
    gtg = (P' * P) .* (B' * B); % compute G'G

    G = S * gtg - mttkrp(X, {B, P, S}, 3) ...
        + alpha * L * S + beta * D' * (D * S * W - Y) * W';

    F = norm(S * khatrirao(B, P)' - ...
        reshape(X.data, size(X, 3), []), 'fro') ^ 2 + ...
        alpha * trace(S' * L * S) + ...
        beta * norm(D * S * W - Y, 'fro') ^ 2;
end