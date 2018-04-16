function [D]=mahal_interdist(X,Alpha)
%MAHAL_ROBUST Robust Mahalanobis distance squared.
% X - Control data matrix (row: samples, column: variables)
% Y - Case data matrix (row: samples, column: variables)
if nargin<2, Alpha=0.75; end
[sig] = robustcov(X,'OutlierFraction',1-Alpha);
n=size(X,1);
D=zeros(n);
for i=1:n-1
  for j=i+1:n
      D(i,j)=(X(i,:)-X(j,:))*(sig\(X(i,:)-X(j,:))');
  end
end
D=D+D';
end
