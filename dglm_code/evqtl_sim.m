N=10000;
testn=200:50:1000;  %sample size
C=zeros(size(testn));

for t=1:length(testn)
    c=0;
for k=1:100
    naa=1600; nAa=4800; nAA=3600;
    y1=[1*randn(naa,1); 2+sqrt(2)*randn(nAa,1); 4+sqrt(3)*randn(nAA,1)];
    x=[zeros(naa,1); ones(nAa,1); 2*ones(nAA,1)];

    n=testn(t);
    idx=randperm(N);
    i=idx(1:n);
    y1=y1(i);
    x=x(i);

    %boxplot(y1,x) % evqtlplot(y1,x,1)
    %{
    mdl = LinearModel.fit(x,y1);
    y2=mdl.Residuals.Raw;
    mdl2 = LinearModel.fit(x,y2.^2);
    mdl2.Coefficients.Estimate(2)    % snp_eff_dispertion
    mdl2.Coefficients.pValue(2)      % snp_dispertion p-value =0.8328
    %}

    %%
    [b]=regress(y1,[ones(size(x)) x]);
    y2=y1-(b(1)+b(2)*x);
    [effs, ~, ~, ~, Stats]=regress(y2.^2,[ones(size(x)) x]);
    p=Stats(3);
    %effsize=effs(2)
    if p<=0.05/n, c=c+1; end
   
end
    C(t)=c/100;
end

%%
h=figure;
plot(testn,C,'o');
hold on
p=polyfit(testn,C,3);
xt=min(testn):max(testn);
f=polyval(p,xt);
plot(xt,f,'r-');
%ylim([0 1])
pubgraph(h)

