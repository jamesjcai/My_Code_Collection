% X=[geno012 covariates];
% [~,~,r]=regress(y1,[ones(size(X,1),1) X]);
% [effs, ~, ~, ~, Stats]=regress(r.^2,[ones(size(geno012,1),1) geno012]);
% p=Stats(3)
% effsize=effs(2)

i=~isnan(x);
x=x(i);
y1=y1(i);
x1=[ones(size(x,1),1) x];

%% method 1 ---
[b,~,y2]=regress(y1,x1);
[bls, ~, ~, ~, Stats]=regress(y2.^2,x1);
p=Stats(3)           % snp_dispertion p-value =0.8328 = mdl2.Coefficients.pValue(2)  
effsize=bls(2)    % snp_eff_dispertion = mdl2.Coefficients.Estimate(2) 

%% method 2
% [b,stats] = robustfit(x1,y1,'cauchy')
[b,~,y2]=regress(y1,x1);
[brob,Stats] = robustfit(x,y2.^2,'cauchy');
% By default, robustfit adds a first column of 1s to X. 
p=Stats.p(2)
effsize=brob(2)


scatter(x,y2.^2,'filled'); grid on; hold on
plot(x,bls(1)+bls(2)*x,'r','LineWidth',2);
plot(x,brob(1)+brob(2)*x,'g','LineWidth',2)
legend('Data','Ordinary Least Squares','Robust Regression')

