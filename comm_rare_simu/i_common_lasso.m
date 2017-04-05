addpath('../mymatlabcentral/');

expr=y;
geno=g012;

figure;
evqtlplot(expr,geno(:,idx_best));

%% 
disp('Using GLMNET method...');
addpath('../PrediXcan_DBPipeline_elasticNet/glmnet/glmnet_matlab/');

family='deviance';
options=[];
type='deviance';
nfolds=10;
foldid=[];
parallel=false;
keep=true;

cvfit = cvglmnet(geno, expr,'gaussian', options, type, nfolds, foldid,...
        parallel, keep);
    
cvglmnetPlot(cvfit)
%%
% Pull info from fit to find the best lambda 
fitdf = [cvfit.cvm cvfit.lambda];
[cvmbest,nrowbest]=min(cvfit.cvm);
bestlam=fitdf(nrowbest,:);

ret=cvfit.glmnet_fit.beta(:,nrowbest);

if any(ret~=0)
    fprintf('Good - predictable\n');
    predmat=cvfit.fit_preval(:,nrowbest);
    mdl = fitlm(predmat,expr);
    rsq=mdl.Rsquared.Ordinary
    pval=mdl.Coefficients.pValue(2)
    % f = regstats(expr,predmat,'linear','fstat');
    % pval2=f.fstat.pval
    figure;
        [expr,i]=sort(expr);
        predmat=predmat(i);
        c = linspace(1,10,length(predmat));
        scatter(predmat,expr,25,c,'filled')
        xlabel('y predicted by using genotype')
        xlabel('y observed')
        box on
        
     figure;
     bar(ret);
     hold on
     line([idx_best idx_best],ylim,'color','r')
     title('beta estimated using GLMNET method')
else
    fprintf('Bad - unpredictable\n');
end

%% Using naive LASSO method
%{
disp('Using naive LASSO method...');
[B, FitInfo] = lasso(geno,expr,'CV',10);
lassoPlot(B,FitInfo,'PlotType','CV');

sum(B(:,FitInfo.IndexMinMSE)>0)
figure;
bar(B(:,FitInfo.IndexMinMSE))
hold on
line([idx_best idx_best],ylim,'color','r')
title('beta estimated using naive LASSO method')
%}
