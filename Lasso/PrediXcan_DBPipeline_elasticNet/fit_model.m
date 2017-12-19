addpath('../glmnet_matlab');

%geno=csvread('exampledata_unpredictable/cisgenos.txt',1,1);
%expr=csvread('exampledata_unpredictable/exppheno.txt',1,1);

geno=csvread('exampledata_predictable/cisgenos.txt',1,1);
expr=csvread('exampledata_predictable/exppheno.txt',1,1);

%%

% fit=glmnet(geno, expr);
% glmnetPrint(fit)
% glmnetPlot(fit)

%% 
family='deviance';
options=[];
type='deviance';
nfolds=10;
foldid=[];
parallel=false;
keep=true;

cvfit = cvglmnet(geno, expr,'gaussian', options, type, nfolds, foldid,...
        parallel, keep);
    
% cvglmnetPlot(cvfit)
%%
% Pull info from fit to find the best lambda 
fitdf = [cvfit.cvm cvfit.lambda];
[cvmbest,nrowbest]=min(cvfit.cvm);
bestlam=fitdf(nrowbest,:);
ret=cvfit.glmnet_fit.beta(:,nrowbest);
% ret(ret==0)=NaN;

if any(ret~=0) 
    fprintf('Good - predictable\n');
    predmat=cvfit.fit_preval(:,nrowbest);
    mdl = fitlm(predmat,expr);
    rsq=mdl.Rsquared.Ordinary
    pval=mdl.Coefficients.pValue(2)
    % f = regstats(expr,predmat,'linear','fstat');
    % pval2=f.fstat.pval
    scatter(predmat,expr);    
else
    fprintf('Bad - unpredictable\n');
end
