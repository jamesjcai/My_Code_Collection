function [rmse, bt] = l0glmcv(X, y, lambda, cv, epn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-validation with generalized linear model
% Written by Zhenqiu Liu
% Cedars-Sinai Medical Center
%  12/18/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5,
    epn = 1e-4;
end

if nargin < 4,
    cv = 5;
end

[n, m] = size(X);
     

cvidx = crossvalind('Kfold', n, cv);
y = y(:);
%cvidx = crossvalid('Kfold', y, cv); % when y is binary
mx  = max(cvidx);
Re = [];
for i =1:mx,
     te = (cvidx ==i);
     tr = ~te;
     Xr = X(tr,:);
     yr = y(tr);
     Xt = X(te,:);
     yt = y(te);
     w = lognewtonl0f(Xr, yr, lambda);
     w(abs(w) < epn) = 0;
     E = logistestt(Xt, yt, w);
     Re = [Re; E];
end
rmse = sqrt(mean(Re.^2));
bt = lognewtonl0f(X, y, lambda);
bt(abs(bt) < epn) = 0;
end

