function [rd]=mahal_robust(Y,X,alphav,weight,methodid)
%MAHAL_ROBUST Robust Mahalanobis distance squared.
% X - Control data matrix (row: samples, column: variables)
% Y - Case data matrix (row: samples, column: variables)
if nargin<5, methodid=1; end
if nargin<4, weight=[]; end
if nargin<3 || isempty(alphav), alphav=0.75; end
switch methodid
    case 1
        % new function since v2016a. 2x faster than mcdcov
        % requirs: size(X,2)<2*size(X,1)
        [sig,mu] = robustcov(X,'OutlierFraction',1-alphav);
    case 2
        % LIBRA toolbox requirs: size(X,2)<=50 
        rew=mcdcov(X,'plots',0,'alpha',alphav);
        mu=rew.center;
        sig=rew.cov;
end
% for i=1:n
%      % rd(i)=(Y(i,:)-mu)*inv(sig)*(Y(i,:)-mu)';
%      rd(i)=(Y(i,:)-mu)*(sig\(Y(i,:)-mu)');
% end
if ~isempty(weight)
    S=weight.*pinv(sig).*weight;
else
    S=pinv(sig);
end
rd=(Y-mu)*S*(Y-mu)';
rd=diag(rd);
