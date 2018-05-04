
rng('default');
load('C:\Users\jcai\Documents\GitHub\My_Code\levene_tests\testdata_scz_Z1_ctl_Z0.mat')
mu0=mean(Z0);
sigma0=cov(Z0);
mu1=mean(Z1);
sigma1=cov(Z1);

R0=mvnrnd(mu0,sigma0,212);
R1=mvnrnd(mu1,sigma1,214);

[p]=anderson_2006_test(R0,R1,0.5,'on')



