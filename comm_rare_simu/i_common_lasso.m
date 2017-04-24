

% 
% % addpath('../mymatlabcentral/');
addpath('supporting/');
figure;
subplot(2,2,1)
evqtlplot(expr,g012(:,idx1));
xlabel(D{1}{1})

subplot(2,2,2)
evqtlplot(expr,g012(:,idx2));
xlabel(D{1}{2})

% if length(D{1})>2
%     subplot(2,2,3)
%     evqtlplot(expr,g012(:,idx3));
%     xlabel(D{1}{3})
% end
%% 
disp('Using GLMNET method...');
addpath('../glmnet_win64_Selmaan');
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

idx_comm=find(maf>=0.05);
geno_comm=snp_pickmarker(geno,[],idx_comm);
g012_comm=snp_012geno(geno_comm);


cvfit = cvglmnet(g012_comm, expr,'gaussian', options, typet, nfolds, foldid,...
        useparallel, keep);

[xl,yl,xs,ys,beta,pctvar,mse] = plsregress(g012_comm,expr,4);    
beta=beta(2:end);

    
% cvglmnetPlot(cvfit)
%%
% Pull info from fit to find the best lambda 
fitdf = [cvfit.cvm cvfit.lambda];
[cvmbest,nrowbest]=min(cvfit.cvm);
bestlam=fitdf(nrowbest,:);

ret=cvfit.glmnet_fit.beta(:,nrowbest);

if any(ret~=0)
    fprintf('Good - predictable\n');
    predmat=cvfit.fit_preval(:,nrowbest);
    
    mdl=fitlm(predmat,expr);
    rsq=mdl.Rsquared.Ordinary
    pval=mdl.Coefficients.pValue(2);
    
    % f = regstats(expr,predmat,'linear','fstat');
    % pval2=f.fstat.pval
    %figure;
    subplot(2,2,4)
        g012_rare=g012(:,idx3);

%        [expr,i]=sort(expr);
%        predmat=predmat(i);

%         c = linspace(1,10,length(predmat));
%         scatter(predmat,expr,25,c,'filled')
        hold on
        scatter(predmat,expr);
        lsline;
%        plot(predmat(g012_rare==1),expr(g012_rare==1),'x','markersize',15)
%        plot(predmat(g012_rare==2),expr(g012_rare==2),'+','markersize',15)
        
        xlabel('y predicted by GLMNET')
        ylabel('y observed')
        box on
        
        expr_residual=expr-predmat;

	
  

             
        
     figure;

     subplot(3,1,1)
     retv=zeros(size(geno,2)/2,1);
     retv(idx_comm)=beta;
     bar(retv);

     subplot(3,1,2)
     retv=zeros(size(geno,2)/2,1);
     retv(idx_comm)=ret;
     bar(retv);
     
     [~,idx_B]=sort(abs(retv),'descend');
     hold on

     line([idx_B(1) idx_B(1)],ylim,'color','r','linestyle','-')
     line([idx_B(2) idx_B(2)],ylim,'color','m','linestyle','-')
     line([idx_B(3) idx_B(3)],ylim,'color','g','linestyle','-')
     
     line([idx_best(1) idx_best(1)],ylim,'color','r','linestyle','--')
     line([idx_best(2) idx_best(2)],ylim,'color','m','linestyle','--')
     line([idx_best(3) idx_best(3)],ylim,'color','g','linestyle','--')
     
     title('Beta estimated using GLMNET')
     % xlabel(['r:',D{1}{1},' m:',D{1}{2},' g:',D{1}{3}])
     
     subplot(3,1,3)
     bar(maf)
     hold on
     line([idx1 idx1],ylim,'color','r')
     line([idx2 idx2],ylim,'color','m')
%     line([idx3 idx3],ylim,'color','g')
     
     line(xlim,[0.05 0.05],'color','r');
     ylim([0 0.5]);
     title('MAF and causal mutations (r,m,g)')
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
