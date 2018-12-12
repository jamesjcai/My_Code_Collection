function [lgu,dropr,lgcv,glist]=i_scmatrix_stat(X,glist)

if nargin<2, glist=[]; end

dropr=1-sum(X>0,2)./size(X,2);
u=nanmean(X,2);
cv=nanstd(X,[],2)./u;
lgu=log10(u);
lgcv=log10(cv);

i=isnan(lgu) | isinf(lgu) | isnan(lgcv) | isinf(lgcv);
lgu(i)=[];
lgcv(i)=[];
dropr(i)=[];
if ~isempty(glist), glist(i)=[]; end

[xyz,i]=sortrows([lgu dropr lgcv],[1 2 3]);
lgu=xyz(:,1);
dropr=xyz(:,2);
lgcv=xyz(:,3);
if ~isempty(glist), glist=glist(i); end
