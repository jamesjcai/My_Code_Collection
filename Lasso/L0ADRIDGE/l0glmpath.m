function [outopt, E, B] = l0glmpath(X, y, lambda,  cv, epn)

if nargin < 5,
    epn = 1e-4;
end

if nargin < 5,
    cv = 4;
end

m = size(X,2);

E = zeros(length(lambda), 1);

B = zeros(m, length(lambda));

for i = 1:length(lambda),
    [rmse, bt] = l0glmcv(X, y, lambda(i), cv, epn);
    B(:, i) = bt;
    E(i) = rmse;
end
 [e, j] = min(E);
 outopt.e = e;
 outopt.b =B(:, j);
 outopt.lam = lambda(j);
end

