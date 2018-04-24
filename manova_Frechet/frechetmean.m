function [u]=frechetmean(X,dist)

if nargin<2, dist = 'euc'; end

switch dist
    case {'euc','euclidean'}
        [~,u]=kmeans(X,1);
    case {'mah','mahal','mahalanobis'}
        disp('xxx')
        %S=nancov(X);        
        [S,u0] = robustcov(X,'OutlierFraction',0.25);
        %u0=nanmean(X);
        %u0=ones(size(u0));
        options = optimset('fminsearch');
        % u=fminsearch(@i_mahal,u0,options,X,S);
        u=fminsearchbnd(@i_mahal,u0,...
                zeros(size(u0)),inf(size(u0)),options,X,S);
        % trace((X-mu)*inv(S)*(X-mu)')
end

end

function [d]=i_mahal(x,X,S)
    x=abs(x);
    d=trace((X-x)*(S\(X-x)'));
end
        
