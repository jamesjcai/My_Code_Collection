load testdata_scz_Z1_ctl_Z0.mat
[~,score]=pca([Z0;Z1]);
figure;
plot(score(1:212,1),score(1:212,2),'b.','markersize',15)
hold on
plot(score(213:end,1),score(213:end,2),'r.','markersize',15);
