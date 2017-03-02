function [rd]=mahal_robust(Y,X,Alpha,Methodid)
% X - Control data matrix (row: samples, column: variables)
% Y - Case data matrix (row: samples, column: variables)
if nargin<4, Methodid=1; end
if nargin<3, Alpha=0.75; end

n=size(Y,1);
rd=nan(n,1);
switch Methodid
    case 1
        % new function since v2016a. 2x faster than mcdcov
        % requirs: size(X,2)<2*size(X,1)
        [sig,mu] = robustcov(X,'OutlierFraction',1-Alpha);
    case 2
        % requirs: size(X,2)<=50
        rew=mcdcov(X,'plots',0,'alpha',Alpha);
        mu=rew.center;
        sig=rew.cov;
end
for i=1:n
     % rd(i)=(Y(i,:)-mu)*inv(sig)*(Y(i,:)-mu)';
     rd(i)=(Y(i,:)-mu)*(sig\(Y(i,:)-mu)');
end
