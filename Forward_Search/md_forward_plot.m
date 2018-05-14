function [M,d_min]=md_forward_plot(X)
%load testdata_scz_Z1_ctl_Z0.mat
%X=Z1;
m=size(X,1);
m2=round(m/2);

[~,~,mah]=robustcov(X,'OutlierFraction',0.5);
[~,idx]=sort(mah);

M=[];
d_min=min(mah(idx(m2+1:end)));
CovMat=[];
for k=1:m-m2
    d_min=[d_min min(mah(idx(m2+k+1:end)))];
    Zm=X(idx(1:m2+k),:);
    CovMat=[CovMat i_triu2vec(cov(Zm))];
    mah=sqrt(mahal(X,Zm));
    [~,idx]=sort(mah);
    M=[M mah];
end
figure;
subplot(2,2,1)
plot(M','k')
xlabel('Subset size m')
ylabel('MD')
subplot(2,2,2)
plot(d_min,'-k');
xlabel('Subset size m')
ylabel('Minimum MD')
subplot(2,2,3)
plot(CovMat','k')
xlabel('Subset size m')
ylabel('Elements of Cov.')
return;
%%
addpath('../FSDA/FSDA');
addpath('../FSDA/FSDA/clustering/');
addpath('../FSDA/FSDA/combinatorial/');
addpath('../FSDA/FSDA/graphics/');
addpath('../FSDA/FSDA/multivariate/');
addpath('../FSDA/FSDA/regression/');
addpath('../FSDA/FSDA/utilities/');
addpath('../FSDA/FSDA/utilities_stat/');
[out]=FSMeda(X,0);
fground=struct;
fground.flabstep='';
malfwdplot(out,'fground',fground);







