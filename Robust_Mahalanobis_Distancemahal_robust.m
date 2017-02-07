function [rd]=mahal_robust(X,Y,methodid)
% control expr: X; case expr: Y (row: sample)
if nargin<3, methodid=1; end
n=size(Y,1);
rd=nan(n,1);    
switch methodid
    case 1
        % new function since v2016a. 2x faster than mcdcov
        [sig,mu] = robustcov(X);        
    case 2
        rew=mcdcov(X,'plots',0);
        mu=rew.center;
        sig=rew.cov;
end
for i=1:n
    rd(i)=sqrt((Y(i,:)-mu)*inv(sig)*(Y(i,:)-mu)');
end

