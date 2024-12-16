function [N,p,X,y] = iSyntheticData1
rng default
N = 10000;
p = 30;
useful = [6,8,9,11,13];
C = randn(p,p);
R = corrcov(C'*C);
X = mvnrnd(zeros(p,1),R,N);
% Make features 15 to 19 highly correlated with useful features:
% 15 -> 6
% 16 -> 8
% 17 -> 9
% 18 -> 11
% 19 -> 13
corrStd = 0.1;
X(:,15:19) = X(:,useful) + corrStd*randn(N,5);
noiseStd = 0.1;
t = 0.5*cos(X(:,11)) + sin(X(:,9).*X(:,8)) + 0.5*X(:,13).*X(:,6) + noiseStd*randn(N,1);
y = rescale(t,0,1);
X = zscore(X);
end