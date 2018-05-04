rng('default');
load('C:\Users\jcai\Documents\GitHub\My_Code\levene_tests\testdata_scz_Z1_ctl_Z0.mat')
mu0=mean(Z0);
sigma0=cov(Z0);
mu1=mean(Z1);
sigma1=cov(Z1);
R0=mvnrnd(mu0,sigma0,212);
R1=mvnrnd(mu1,sigma1,214);
%%
R=[R0; R1];
grp=[zeros(212,1);ones(214,1)];

[B W]=scattermat(R,grp)

%%
N=[212 214];
k=2;
x=R';
w = 0;
T=1;
for i = 1:k
    for j = 1:N(i)
        x_diff = x(i,j) - mean(x(i,:));
        w = w + (x_diff * x_diff^T);
    end
end

b = 0;

for i = 1:k
    x_diff = mean(x(i,:)) - mean(x);
    b = b + (N(i) * x_diff * x_diff^T);
end








