function [v]=pv(X)
% Proportional Variability
%
% http://onlinelibrary.wiley.com/doi/10.1111/j.2006.0030-1299.15067.x/epdf
% http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0084074
% x=[100*ones(25,1); 10*ones(50,1); 100*ones(25,1)];
% n=numel(x);

%dv=abs(bsxfun(@minus,x,x'));
%dv2=bsxfun(@max,x,x');
% pv=2*sum(sum(abs(x-x')./max(x,x')))./(numel(x)*(numel(x)-1));

[m,n]=size(X);
v=zeros(m,1);

for k=1:m
    x=X(k,:);
    d=1-min(x,x')./max(x,x');
    %d(isnan(d))=0;
    v(k)=nansum(d(:))./(n*(n-1));
end



