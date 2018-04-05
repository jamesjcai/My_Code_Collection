function [s]=pca_sim_factor(Xs,Xh)

L = pca(Xs);
M = pca(Xh);

k=3;
L=L(:,1:k);
M=M(:,1:k);

s=trace(L'*M*(M')*L)/k;

