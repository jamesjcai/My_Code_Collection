% Source code by Gavin Lloyd and Richard Brebreton.

% create matrix of random nos
% J >> I
X=randn(10,1000);
% mean centre
X=bsxfun(@minus,X,mean(X));

%% standard mahalanobis distance
C=cov(X);
IC=inv(C); % should generate 'matrix close to singular' warning
for i=1:size(X,1)
   MS(i,1)=X(i,:)*IC*X(i,:)'; 
end

%% use PCA
[U,S,V]=svd(X); % this is equivalent to PCA
T=U*S; % PCA scores. V are loadings. diag(S)^2 = eigenvalues
% NB only min(I,J) components will have eigenvalue>0
lo=min([size(X,1),size(X,2)]);
T=T(:,1:lo); % or you could reduce to a smaller number of components if you wanted

% autoscale scores
M=mean(T);
S=std(T);
Ts=bsxfun(@minus,T,M); % NB scores already centred if X is centred
Ts=bsxfun(@rdivide,Ts,S);

% covariance of T
CT=cov(T); % should be identity matrix because PCA makes components orthogonal i.e. no covariance between components

% mahal dist
% since CT is identity matrix we can ignore it in the mahal calc (equiv divide by 1)
for i=1:size(X,1)
    MD(i,:)=T(i,:)*T(i,:)';
end
