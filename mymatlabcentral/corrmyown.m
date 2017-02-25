function [r,p]=corrmyown(a,b,c,quietpls,plotit)
if nargin<5
    plotit=false;
end
if nargin<4
    quietpls=false;
end
if nargin<3
    c=[];
end
if size(a,2)~=1&&size(a,1)==1, a=a'; end
if size(b,2)~=1&&size(b,1)==1, b=b'; end


if isempty(c)
    idx=~(isnan(a)|isnan(b));
    [r,p]=corr(a(idx),b(idx),'type','s');
    [r1,p1]=corr(a(idx),b(idx),'type','p');
else
    idx=~(isnan(a)|isnan(b)|isnan(c));
    [r,p]=partialcorr(a(idx),b(idx),c(idx),'type','s');
    [r1,p1]=partialcorr(a(idx),b(idx),c(idx),'type','p');

fprintf('*Partial correlations:\n');
end
    if ~quietpls
        fprintf('Spearman''s rho = %.3f, P = %g (n = %d)\n',r,p,length(a(idx)));
        fprintf('Pearson''s r = %.3f, P = %g\n',r1,p1);
        fprintf('\n%.3f (%g)\t',r,p);
        fprintf('%.3f (%g)\n\n',r1,p1);
    end
    if plotit
        scattercloud(a(idx),b(idx));
        addreg(a(idx),b(idx));
        box on
        title(sprintf('Spearman''s rho = %.3f, P = %g\nPearson''s r = %.3f, P = %g',...
              r,p,r1,p1))        
    end
