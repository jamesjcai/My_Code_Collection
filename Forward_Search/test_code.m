load testdata_scz_Z1_ctl_Z0.mat


m=size(Z0,1);
M=[];
for k=0.5:-0.01:0
    [~,~,mah1] = robustcov(Z0,'OutlierFraction',k);
    M=[M mah1];
end
figure; 
plot(M');
