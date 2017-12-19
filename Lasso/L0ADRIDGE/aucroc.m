function [auc, tp, fp] = aucroc(y, p)
% This function calculated the area under curve
% y: the true labels
% p: the probability prediction or (score function)
% auc: the area under curve
% tp: the true positive rate
% fp: the false positive rate
% 7/15/2006
% Zhenqiu Liu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% process targets

y = y > 0;   % transfer the label to 0/1 if not.

% sort by classifier output

[P,idx] = sort(-p);
y       = y(idx);

% compute true positive and false positive rates

tp = cumsum(y)/sum(y);
fp = cumsum(~y)/sum(~y);

% add trivial end-points

tp = [0 ; tp ; 1];
fp = [0 ; fp ; 1];

n = size(tp, 1);
auc = sum((fp(2:n) - fp(1:n-1)).*(tp(2:n)+tp(1:n-1)))/2;