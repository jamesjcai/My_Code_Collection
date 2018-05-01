function [p,stats]=anderson_2006_test(Z0,Z1,alphav)
% Distance?Based Tests for Homogeneity of Multivariate Dispersions
% MARTI J. ANDERSON

if nargin<3, alphav=0.75; end
n0=size(Z0,1);
n1=size(Z1,1);
grp=[zeros(n0,1); ones(n1,1)];

% Z0=zscore(Z0);
% Z1=zscore(Z1);
[~,~,mah0] = robustcov(Z0,'OutlierFraction',1-alphav);
[~,~,mah1] = robustcov(Z1,'OutlierFraction',1-alphav);
z=[mah0;mah1];
[p,~,stats]=anova1(z,grp,'off');       
