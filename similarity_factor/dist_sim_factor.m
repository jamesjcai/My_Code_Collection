function [s,d]=dist_sim_factor(Xs,Xh)

u=mean(Xh)-mean(Xs);
[sig]=robustcov(Xs);
%d=sqrt(u*pinv(sig)*u');
d=sqrt(u*(sig\u'));

% [~,S,~] = svd(Xs); % A = U*S*V'
% T=S;
% T(S~=0) = 1./S(S~=0);
% %svdInvA =V*T'*U';
% sig=T(1:min(size(Xs)),1:min(size(Xs)))
% d=sqrt(u*sig*u');

fun = @(z)exp(-z.^2/2);
m=sqrt(2/pi);
s = m*integral(fun,d,Inf);
