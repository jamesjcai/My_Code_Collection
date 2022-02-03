clear all, close all, clc
load ../DATA/FLUIDS/CYLINDER_ALL.mat
Y = VORTALL + .01*randn(size(VORTALL));

[U,S,V] = svd(Y,'econ');     % SVD matrix
beta = size(Y,2)/size(Y,1);  % aspect ratio of matrix
sigma = diag(S);             % extract singular values
% optimal threshold tau
tau = optimal_SVHT_coef(beta,0)*median(sigma);  

semilogy(sigma,'k-o','LineWidth',1.2)
hold on
semilogy(sigma(sigma>tau),'ro','LineWidth',1.2)
axis([0 150 1 10000])
grid on
set(gcf,'Position',[100 100 250 250])
semilogy(15,sigma(15),'b.','MarkerSize',30)
semilogy(27,sigma(27),'b.','MarkerSize',30)
set(gcf,'PaperPositionMode','auto') % 
% print('-depsc2', '-loose', './sigmavals.eps'); % eps are vector images

%%
plotCylinder(reshape(U(:,27),m,n),n,m)



%%
Y = .01*randn(size(VORTALL));
[U,S,V] = svd(Y,'econ');
subplot(1,2,1)
semilogy(diag(S),'b')
hold on
subplot(1,2,2)
plot(cumsum(diag(S)/sum(diag(S))))
hold on
%%
sigs = diag(S);
beta = size(VORTALL,2)/size(VORTALL,1);
thresh = optimal_SVHT_coef(beta,0)*median(sigs);
subplot(1,2,1)
semilogy(sigs(sigs>thresh),'ro','LineWidth',1.2)