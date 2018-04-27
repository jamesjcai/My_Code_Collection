function [p,F]=anderson_2001_test(Z0,Z1,alphav)
% A new method for non-parametric multivariate analysis of variance
% MARTI J. ANDERSON

if nargin<3, alphav=[]; end
n0=size(Z0,1);
n1=size(Z1,1);
a=[-1*ones(n0,1); ones(n1,1)];
M=a*a';
M(M==-1)=0;
% n=sum(triu2vec(M)>0)
n=0.5*n0*(n0-1)+0.5*n1*(n1-1);
Z=[Z0;Z1];
if ~isempty(alphav)
    [D]=mahal_interdist(Z,alphav);
else
    [D]=interdist(Z);       % Euclidean distance
end
SSt=mean(triu2vec(D).^2);
SSw=sum(triu2vec(D.*M).^2)/n;
SSa=SSt-SSw;
a=2;
N=n0+n1;

df1=a-1; df2=N-a;
F=(SSa/df1)/(SSw/df2);
p=1-fcdf(F,df1,df2);

% y=triu2vec(D).^2;
% SSt=mean(y);
% g=findgroups(triu2vec(M));
% SSw=splitapply(@mean,y,g)
% SSw=SSw(2);
