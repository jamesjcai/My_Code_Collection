function [yes]=i_95pct_error_ellipse_OLD(datax,xy,conf)

if nargin<3, conf=0.95; end    

% Calculate the eigenvectors and eigenvalues
covariance = cov(datax);
[~, eigenval] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c] = find(eigenval == max(max(eigenval)));
%largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    %smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    %smallest_eigenvec = eigenvec(1,:);
end
% Get the 95% confidence interval error ellipse
% chisquare_val = 2.4477;
chisquare_val = sqrt(chi2inv(conf,2));

a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);
% yes=((xy(1)/a)^2 + (xy(2)/b)^2) <= 5.991;

yes=((xy(1)/a)^2 + (xy(2)/b)^2) <= chisquare_val^2;
