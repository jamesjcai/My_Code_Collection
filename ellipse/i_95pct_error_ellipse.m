function [yes]=i_95pct_error_ellipse(datax,xy,conf)

if nargin<3, conf=0.95; end    

% Calculate the eigenvectors and eigenvalues
covariance = cov(datax);
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c,r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

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


% Get the coordinates of the data mean
avg = mean(datax);
% Get the 95% confidence interval error ellipse
% chisquare_val = 2.4477;
chisquare_val = sqrt(chi2inv(conf,2));

a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);
% yes=((xy(1)/a)^2 + (xy(2)/b)^2) <= 5.991;
% yes=((xy(1)/a)^2 + (xy(2)/b)^2) <= chisquare_val^2;
%yes=(((xy(1)-mm(1))/a)^2 + ((xy(2)-mm(2))/b)^2) <= 1;

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the 95% confidence interval error ellipse
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
avg1=avg*R;
X01=avg1(1);
Y01=avg1(2);
xynew=xy*R;
yes=(((xynew(1)-X01)/a)^2 + ((xynew(2)-Y01)/b)^2) <= 1;

