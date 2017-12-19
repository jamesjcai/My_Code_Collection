function [w, u] = l0regprimal(X, y, lambda, p, MaxIt, epn)
% This function computes the Least Squares parameters
% with a penalty on the Lp-norm of the parameters
%
% Method used:
% Primal  EM algorithm
%
% Zhenqiu Liu
% Cedars Sinai Medical Center
% 05/31/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6,
     epn = 1e-6;
end

if nargin < 5,
   MaxIt  = 1000;
end

if nargin < 4,
    p = 1;
end

[n, m] = size(X);

w = (X'*X + lambda*eye(m))\(X'*y);

i = 0;
while  i < MaxIt
    i = i + 1;
    u = w;
    up = abs(u').^(2-p);
    Xu = X.*repmat(up, n, 1);
    XXfull = Xu'*X;
    Xyfull = Xu'*y;
    w = (XXfull + lambda*eye(m))\Xyfull;
    %w(w <0) = 0;
    if norm(w-u) < epn, break; end
    if i > MaxIt,
       disp('Warning: No enough iterations');
       break;
   end
end
w(abs(w) < 1e-3) =0;   