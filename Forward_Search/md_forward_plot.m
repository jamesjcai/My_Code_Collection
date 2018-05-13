% function [M]=md_forward_plot(X)
load testdata_scz_Z1_ctl_Z0.mat
X=Z1;
m=size(X,1);
m2=round(m/2);

[~,~,mah]=robustcov(X,'OutlierFraction',0.5);
[~,idx]=sort(mah);

M=[];
for k=1:m-m2
    Zm=X(idx(1:m2+k),:);
    mah=sqrt(mahal(X,Zm));
    [~,idx]=sort(mah);
    M=[M mah];
end
figure;
plot(M')
