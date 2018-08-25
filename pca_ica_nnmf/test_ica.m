addpath('pca_ica')
load testdata_scz_Z1_ctl_Z0.mat
r=3;

%Zmixed=Z0';
Zmixed=Z0;
Zfica0 = fastICA(Zmixed,r);
Zkica0 = kICA(Zmixed,r);
Zpca0 = PCA(Zmixed,r);

%Zmixed=Z1';
Zmixed=Z1;
Zfica1 = fastICA(Zmixed,r);
Zkica1 = kICA(Zmixed,r);
Zpca1 = PCA(Zmixed,r);

%%
figure;
% Fast ICA
subplot(6,1,1);
for i = 1:r
    plot(Zfica0(i,:),'-'); hold on;
end
title('Independent components [Fast ICA]');
axis tight;

% Max-kurtosis
subplot(6,1,2);
for i = 1:r
    plot(Zkica0(i,:),'-'); hold on;
end
title('Independent components [max-kurtosis]');
axis tight;

% PCA
subplot(6,1,3);
for i = 1:r
    plot(Zpca0(i,:),'-'); hold on;
end
title('Principal components');
axis tight;



% Fast ICA
subplot(6,1,4);
for i = 1:r
    plot(Zfica1(i,:),'-'); hold on;
end
title('Independent components [Fast ICA]');
axis tight;

% Max-kurtosis
subplot(6,1,5);
for i = 1:r
    plot(Zkica1(i,:),'-'); hold on;
end
title('Independent components [max-kurtosis]');
axis tight;

% PCA
subplot(6,1,6);
for i = 1:r
    plot(Zpca1(i,:),'-'); hold on;
end
title('Principal components');
axis tight;

