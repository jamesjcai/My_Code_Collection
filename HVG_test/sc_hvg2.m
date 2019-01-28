function [T]=sc_hvg2(X,genelist,sortit,plotit)
% HVGs selection - This method uses the CV2 on normalized count data to 
% analyze biological variation. 
%
% REF: https://www.nature.com/articles/nmeth.2645
% Input X: Brennecke et al. (2013) uses DESeq's method for normalization.
% 
% USAGE:
% >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >> [X]=sc_norm(X,'type','deseq');
% >> [T]=sc_hvg(X,genelist,true,true);

if nargin<4, plotit=false; end
if nargin<3, sortit=true; end

%[X,sf]=norm_deseq(X);
u=nanmean(X,2);
cv2=nanvar(X,0,2)./u.^2;

xi=1./u; yi=cv2; yi=yi(xi>0.1); xi=xi(xi>0.1);
b=glmfit(xi,yi,'gamma','link','identity');
mdl=fitglm(xi,yi,'linear','Distribution','gamma','link','identity');
cv2fitx=predict(mdl,1./u)

cv2fit=glmval(b,1./u,'identity');
% cv2fit=b(2)./u+b(1);
[cv2fitx cv2fit]

df=size(X,2)-1;
fitratio=cv2./cv2fit;
pval=chi2cdf(fitratio*df,df,'upper');

fdr=mafdr(pval);

T=table(genelist,u,cv2,fitratio,pval,fdr);
T.Properties.VariableNames(1)={'genes'};
if sortit
    T=sortrows(T,'fdr');
end

if plotit
    [~,i]=sort(fitratio,'descend');
    xi=u(i); yi=cv2(i); yifit=cv2fit(i);    
    
    scatter(log(xi),log(yi))
    hold on
    scatter(log(xi(1:10)),log(yi(1:10)),'x');
    plot(log(xi),log(yifit),'.','markersize',10);   
    % plot(log(xi),log(yifit*chi2inv(0.975,df)./df),'.k');
    % plot(log(xi),log(yifit*chi2inv(0.025,df)./df),'.k');
    xlabel('Mean expression, log')
    ylabel('CV^2, log')
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelist(i)};
    end    
    hold off    
end

% Highly variable genes (HVG) is based on the assumption that genes with 
% high variance relative to their mean expression are due to biological 
% effects rather than just technical noise. The method seeks to identify 
% genes that have a higher variability than expected by considering the 
% relationship between variance and mean expression. This relationship is 
% difficult to fit, and in practice genes are ranked by their distance 
% from a moving median (Kolodziejczyk et al., 2015) or another statistic 
% derived from variance is used, e.g. the squared coefficient of variation
% (Brennecke et al. (2013)).
