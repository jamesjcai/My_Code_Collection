function [retval] = wmahalanobis (x, center, cov, weight)
x = (x-mean(x));
%cov = weight * inv(cov);
cov = weight/cov;
retval = diag(x*cov*x');
end