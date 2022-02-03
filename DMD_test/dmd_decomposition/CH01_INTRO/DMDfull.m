function [Phi, omega, lambda, b, Xdmd] = DMDfull(X, varin)
% function [Phi,omega,lambda,b,Xdmd] = DMD(X,varin)
% Computes the Dynamic Mode Decomposition of data X
%
% INPUTS: 
% Columns of X are state snapshots 
% Rows of X are measurements
% Columns X are time points, sampled at equal dt's
%
% Optional parameters: {'parameter_name', [default_value]}
%   {'dt', [1]}
%   {'r', [1e32]}        truncate to rank-r
%   {'nstacks', 1}       number of stacks of the raw data
%
%
% OUTPUTS:
% Phi, the DMD modes
% omega, the continuous-time DMD eigenvalues
% lambda, the discrete-time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstrcted by Phi, omega, b

%% input parsing
p = inputParser; 

% required inputs
p.addRequired('X', @isnumeric);

% parameter value iputs
p.addParameter('dt', 1, @(x)isnumeric(x) && x>0);
p.addParameter('r', 1e32, @(x)isnumeric(x) && x>0);
p.addParameter('nstacks', 1, @(x)isnumeric(x) && x>0);

% now parse the inputs
p.parse(X, varin{:});
inputs = p.Results;

%% stacking the data matrix 
if inputs.nstacks > 1
    Xaug = [];
    for st = 1:inputs.nstacks
        Xaug = [Xaug; X(:, st:end-inputs.nstacks+st)];
    end
    X1 = Xaug(:, 1:end-1);
    X2 = Xaug(:, 2:end);
else
    X1 = X(:, 1:end-1);
    X2 = X(:, 2:end);
end
m = size(X1, 2);

%% DMD
[U, S, V] = svd(X1, 'econ');
r = inputs.r;
if r >= size(U, 2)
    r = size(U, 2);
end

if r >= size(U,2) % no rank truncation
    Atilde = U' * X2 * V / S;
    [W, D] = eig(Atilde); 
    Phi = X2 * V / S * W;
else % truncate modes
    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    Atilde = U_r' * X2 * V_r / S_r;
    [W_r, D] = eig(Atilde);
    Phi = X2 * V_r / S_r * W_r;
end

lambda = diag(D);
omega = log(lambda)/inputs.dt;

%% compute DMD mode amplitudes
x1 = X1(:, 1);
b = Phi\x1;

%% DMD reconstructions
time_dynamics = zeros(r, m);
t = (1:m)/inputs.dt;
for iter = 1:m
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd = Phi * time_dynamics;
