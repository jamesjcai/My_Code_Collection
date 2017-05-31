if ispc
    addpath('../glmnet_win64_Selmaan');
else
    addpath('../glmnet_matlab/');
end

%% 
%{
disp('Using GLMNET method...common snps only');
% addpath('../glmnet_win64');
% addpath('../glmnet_matlab');
% addpath('../PrediXcan_DBPipeline_elasticNet/glmnet/glmnet_matlab/');
% addpath('C:\Users\jcai\Desktop\glmnet\glmnet_matlab');
family='deviance';
options=[];
typet='deviance';
nfolds=10;
foldid=[];
useparallel=false;
keep=true;

idx_comm=find(mafv>=0.05);
geno_comm=snp_pickmarker(geno,[],idx_comm);
g012_comm=snp_012geno(geno_comm);

cvfit = cvglmnet(g012_comm, expr,'gaussian', options, typet, nfolds, foldid,...
        useparallel, keep);
% cvglmnetPlot(cvfit)


%disp('Using PLS method...');    
%[xl,yl,xs,ys,beta,pctvar,mse] = plsregress(g012_comm,expr,4);    
%beta=beta(2:end);
    


% Pull info from fit to find the best lambda 
fitdf = [cvfit.cvm cvfit.lambda];
[cvmbest,nrowbest]=min(cvfit.cvm);
bestlam=fitdf(nrowbest,:);

ret=cvfit.glmnet_fit.beta(:,nrowbest);

if any(ret~=0)
    fprintf('Good - predictable\n');
    predmat=cvfit.fit_preval(:,nrowbest);
    
    mdl=fitlm(predmat,expr);
    rsq=mdl.Rsquared.Ordinary;
    pval=mdl.Coefficients.pValue(2);
    
    % f = regstats(expr,predmat,'linear','fstat');
    % pval2=f.fstat.pval
    figure;
    subplot(2,2,4)
        g012_rare=g012(:,idx3);
%        [expr,i]=sort(expr);
%        predmat=predmat(i);

%         c = linspace(1,10,length(predmat));
%         scatter(predmat,expr,25,c,'filled')
        hold on
        scatter(predmat,expr);
        lsline;
%       plot(predmat(g012_rare==1),expr(g012_rare==1),'x','markersize',15)
%       plot(predmat(g012_rare==2),expr(g012_rare==2),'+','markersize',15)
        xlabel('y predicted by GLMNET')
        ylabel('y observed')
        box on
        expr_residual=expr-predmat;
       
     figure;
    % subplot(3,1,1)
    % retv=zeros(size(geno,2)/2,1);
    % retv(idx_comm)=beta;
    % bar(retv);
     subplot(3,1,1)
     retv=zeros(size(geno,2)/2,1);
     retv(idx_comm)=ret;
     bar(retv);
     
%     [~,idx_B]=sort(abs(retv),'descend');
%     hold on
%     line([idx_B(1) idx_B(1)],ylim,'color','r','linestyle','-')
%     line([idx_B(2) idx_B(2)],ylim,'color','m','linestyle','-')
%     line([idx_B(3) idx_B(3)],ylim,'color','g','linestyle','-')
%     line([idx_best(1) idx_best(1)],ylim,'color','r','linestyle','--')
%     line([idx_best(2) idx_best(2)],ylim,'color','m','linestyle','--')
%     line([idx_best(3) idx_best(3)],ylim,'color','g','linestyle','--')
     title('Beta estimated using GLMNET')
     % xlabel(['r:',D{1}{1},' m:',D{1}{2},' g:',D{1}{3}])
else
    fprintf('Bad - unpredictable\n');
end
%}
%%

disp('Using GLMNET method... all SNPs');
family='deviance';
options=[];
typet='deviance';
nfolds=10;
foldid=[];
useparallel=false;
keep=true;
tic;
cvfit = cvglmnet(g012, expr,'gaussian', options, typet, nfolds, foldid,...
        useparallel, keep);
% Pull info from fit to find the best lambda 
fitdf = [cvfit.cvm cvfit.lambda];
[cvmbest,nrowbest]=min(cvfit.cvm);
bestlam=fitdf(nrowbest,:);
ret=cvfit.glmnet_fit.beta(:,nrowbest);
toc;
if any(ret~=0)
    % fprintf('Good - predictable\n');
    predmat=cvfit.fit_preval(:,nrowbest);
    
    mdl=fitlm(predmat,expr);
    rsq=mdl.Rsquared.Ordinary
    pval=mdl.Coefficients.pValue(2);
    
    % f = regstats(expr,predmat,'linear','fstat');
    % pval2=f.fstat.pval
    figure(fc);
    subplot(2,2,4)
        hold on
        scatter(predmat,expr);
        lsline;
        xlabel('y predicted by GLMNET')
        ylabel('y observed')
        box on
        expr_residual=expr-predmat;
       
     figure(fb);
    % subplot(3,1,1)
    % retv=zeros(size(geno,2)/2,1);
    % retv(idx_comm)=beta;
    % bar(retv);
     subplot(6,1,3)
     bar(ret);
     title('Beta estimated using GLMNET')
     ylabel('Beta')
else
    fprintf('Bad - unpredictable\n');
end

%%
tic;
disp('Using PLS method...');    
[xl,yl,xs,ys,betapls,pctvar,mse] = plsregress(g012,expr,2);    
betapls=betapls(2:end);
figure(fb);
subplot(6,1,4)
bar(betapls)
title('Beta estimated using PLS')
ylabel('Beta')
toc;   

%% Using LASSO method
tic;
disp('Using LASSO method...');
[Bb, FitInfo]=lassoglm(g012,expr,'normal','CV',10);
minpts = find(Bb(:,FitInfo.IndexMinDeviance));
min1pts = find(Bb(:,FitInfo.Index1SE));
%lassoPlot(Bb,FitInfo,'plottype','CV');
figure(fb);
subplot(6,1,5)
bar(Bb(:,FitInfo.IndexMinDeviance))
title('Beta estimated using LASSO')
ylabel('Beta')
toc;

%% Using naive LASSO method
tic;
disp('Using naive LASSO method...');
[B_ln, FitInfo] = lasso(g012,expr,'CV',10);
%lassoPlot(B,FitInfo,'PlotType','CV');
%sum(B(:,FitInfo.IndexMinMSE)>0)
figure(fb);
subplot(6,1,6)
bar(B_ln(:,FitInfo.IndexMinMSE))
title('Beta estimated using naive LASSO')
ylabel('Beta')
toc;