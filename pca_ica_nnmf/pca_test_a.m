X=[2.5 0.5 2.2 1.9 3.1 2.3 2.0 1.0 1.5 1.1; 2.4 0.7 2.9 2.2 3.0 2.7 1.6 1.1 1.6 0.9]';
Xc=X-mean(X);
C=cov(X);
[v,d]=eig(C);   % d = eigenvalues
totalvarexplained = diag(d)./sum(diag(d));

[coeff,score,latent]=pca(X);
coeff(:,1) == v(:,2)
max(coeff(:,1)-v(:,2))

[Xc*coeff Xc/coeff' score]

figure;
subplot(2,2,1)
scatter(X(:,1), X(:,2))
subplot(2,2,2)
X2=(coeff*X')';
scatter(X2(:,1), X2(:,2))

subplot(2,2,3)
X2=X*coeff;
scatter(X2(:,1), X2(:,2))

subplot(2,2,4)
scatter(score(:,1), score(:,2))
% https://www.youtube.com/watch?v=kn_rLM1ZS2Q

dot(coeff(:,1),coeff(:,2))
% https://math.stackexchange.com/questions/403321/which-vectors-are-perpendicular-to-each-other



%%

nm=mean(X,1);  % mean over examples
X0=X-repmat(nm,[n,1]); % subtract mean
S=X0'*X0/n;  % S=cov(X);
[V,D]=eig(S);   % can be slow!
v=V(:,end-k+1:end);  % k largest eigenvectors,, can also fund incrementally ("eigs")


