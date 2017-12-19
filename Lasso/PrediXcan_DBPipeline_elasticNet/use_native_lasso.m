geno=csvread('exampledata_predictable/cisgenos.txt',1,1);
expr=csvread('exampledata_predictable/exppheno.txt',1,1);

[B, FitInfo] = lasso(geno,expr,'CV',10);
lassoPlot(B,FitInfo,'PlotType','CV');

sum(B(:,FitInfo.IndexMinMSE)>0)

figure;
k=FitInfo.IndexMinMSE;
scatter(geno*B(:,k)-FitInfo.Intercept(k),expr)

for k=1:5
    figure; scatter(geno*B(:,10*k),expr)
end
