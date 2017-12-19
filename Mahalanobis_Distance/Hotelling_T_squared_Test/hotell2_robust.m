function [pval,T2,F] = hotell2_robust(X,Y,Alpha)

if nargin<3
    Alpha=0.75;
end

[nx,px] = size(X);
[ny,py] = size(Y);

if px ~= py
   error('# of columns in X and Y must match');
else
   p = px;
end

n = nx + ny;
%mux = mean(x);
%muy = mean(y);

%Sx = cov(x);
%Sy = cov(y);
[Sx,mux] = robustcov(X,'OutlierFraction',1-Alpha);
[Sy,muy] = robustcov(Y,'OutlierFraction',1-Alpha);

% Hotelling T2 statistic, Section 3.6.1 Mardia et al.
    % Su = (nx*Sx + ny*Sy) / (n-2); % unbiased estimate
    Su = ((nx-1)*Sx + (ny-1)*Sy) / (n-2);

d = mux - muy;
D2 = d*inv(Su)*d';
T2 = ((nx*ny)/n)*D2;
F = T2 * (n-p-1) / ((n-2)*p);

pval = 1 - fcdf(F,p,n-p-1);

if nargout == 0
   fprintf('-------------------------------\n');
   fprintf('  nx = %g\n',nx);
   fprintf('  ny = %g\n',ny);
   fprintf('  mean(x) = ');
   fprintf('%1.3f, ',mux);
   fprintf('\n');
   fprintf('  mean(y) = ');
   fprintf('%1.3f, ',muy);
   fprintf('\n');
   fprintf('  T2 = %5.3f\n',T2);
   fprintf('  F(%g,%g) = %5.3f\n',p,n-p-1,F);
   fprintf('  p = %5.5f\n',pval);
   fprintf('-------------------------------\n');
end