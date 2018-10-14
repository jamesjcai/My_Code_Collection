function [cv,robustcv]=mycv(x)
cv=nanstd(x,0,2)./nanmean(x,2);
if nargout>1
    robustcv=iqr(x,2)./nanmedian(x,2);
end
