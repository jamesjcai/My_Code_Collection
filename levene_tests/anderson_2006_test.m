function [p,stats,d,grp]=anderson_2006_test(Z0,Z1,alphav,displayopt)
% Distance-Based Tests for Homogeneity of Multivariate Dispersions
% MARTI J. ANDERSON

if nargin<3 || isempty(alphav), alphav=0.75; end
if nargin<4, displayopt = 'off'; end
n0=size(Z0,1);
n1=size(Z1,1);
grp=[zeros(n0,1); ones(n1,1)];
% Z0=zscore(Z0);
% Z1=zscore(Z1);
%[~,~,mah0] = robustcov(Z0,'OutlierFraction',1-alphav);
%[~,~,mah1] = robustcov(Z1,'OutlierFraction',1-alphav);
[mah0]=sqrt(mahal_robust(Z0,Z0,alphav));
[mah1]=sqrt(mahal_robust(Z1,Z1,alphav));
d=[mah0;mah1];
[p,~,stats]=anova1(d,grp,'off');

if strcmp(displayopt,'on')
    % h=boxplot(d,grp,'notch','on','color','k','symbol','');
    h=boxplot(d,grp,'notch','on','color','k');
    set(h(6,:),'color','k','linewidth',3);
    hold on
    plot(grp+1.25,d,'ok');
end
