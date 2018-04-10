function [s,d]=dist_sim_factor(Xh,Xs,isrobust)
% Xh - historical data (reference/control)
% Xs - snapshot data (test)
% isrobust - false (default) using Fast-MCD to compute the robust covariance matrix

if nargin<3, isrobust=false; end

u=mean(Xh)-mean(Xs);
if isrobust
    sig=robustcov(Xh);    
else
    sig=cov(Xh);
end
d=sqrt(u*(sig\u'));  %d=sqrt(u*pinv(sig)*u');


% [~,S,~] = svd(Xs); % A = U*S*V'
% T=S;
% T(S~=0) = 1./S(S~=0);
% %svdInvA =V*T'*U';
% sig=T(1:min(size(Xs)),1:min(size(Xs)))
% d=sqrt(u*sig*u');

fun = @(z)exp(-z.^2/2);
m=sqrt(2/pi);
s = m*integral(fun,d,Inf);
