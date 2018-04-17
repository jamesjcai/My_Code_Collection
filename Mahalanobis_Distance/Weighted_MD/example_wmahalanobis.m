load('fisheriris.mat')
x = meas(1:50,:);
center = mean(x);
covar = cov(x);
weight = diag(repmat(0.25,1,size(x,2)));
wmahalanobis(x, center, covar, weight)