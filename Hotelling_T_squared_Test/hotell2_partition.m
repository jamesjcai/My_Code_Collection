function [f,c]=hotell2_partition(X,Y)
% http://users.mct.open.ac.uk/paul.garthwaite/HotTwoSamples.html
% The contributions of each variable to Hotelling's T2 statistic. The sum of these contributions equals T2.
%
% F - The contributions of each variable expressed as percentages.
% C - The correlation between each variable and its surrogate: the partition constructs an orthogonal set of surrogate variables that are used to evaluate the individual contributions of variables. Low correlations can be indicative of collinearities between some variables.
%
% There are two independent samples from populations that have multivariate normal distributions. It is assumed that these distributions have the same variance but their means may differ. Hotelling's two-sample T2 is used to test the hypothesis that the two means are equal. The program performs the hypothesis test and determines the contributions of individual variables to Hotelling's T2 statistic.

[nx,px] = size(X);
[ny,py] = size(Y);

if px ~= py
   error('# of columns in X and Y must match');
end

n = nx + ny;
mux = mean(X);
muy = mean(Y);
Sx = cov(X);
Sy = cov(Y);
%[Sx,mux] = robustcov(X,'OutlierFraction',1-Alpha);
%[Sy,muy] = robustcov(Y,'OutlierFraction',1-Alpha);

% Hotelling T2 statistic, Section 3.6.1 Mardia et al.
% Su = (nx*Sx + ny*Sy) / (n-2); % unbiased estimate
Su = ((nx-1)*Sx + (ny-1)*Sy) / (n-2);   % poolvar
df = mux - muy;

% Start that part of the Garthwaite-Koch partition that is the same for all items 
Dmat = diag(1./sqrt(diag(Su)));
DSD=Dmat*Su*Dmat;

[eigvec,eigval]=eig(DSD);
InvRootEig=1./sqrt(diag(eigval));

InvDSDhalf=eigvec*diag(InvRootEig)*eigvec';

W=InvDSDhalf*Dmat*df';
Wsquared = diag(W*W');
GKcontrib =((nx*ny)/n) * Wsquared;
f=GKcontrib./sum(GKcontrib);
c=diag(inv(InvDSDhalf));


