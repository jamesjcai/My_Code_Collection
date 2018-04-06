function [s]=pca_sim_factor(Xs,Xh)

%[coeff1,score1,latent1,tsquared1,explained1] = pca(X1);
[L,~,~,~,explainedL] = pca(Xs);
[M,~,~,~,explainedM] = pca(Xh);

k1=find(cumsum(explainedL)>95,1);
k2=find(cumsum(explainedM)>95,1);

k=max([k1 k2]);
L=L(:,1:k);
M=M(:,1:k);

s=trace(L'*M*(M')*L)/k;

