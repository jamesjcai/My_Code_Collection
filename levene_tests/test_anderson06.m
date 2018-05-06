rng('default');
load('testdata_scz_Z1_ctl_Z0.mat')
mu0=mean(Z0);
sigma0=cov(Z0);
mu1=mean(Z1);
sigma1=cov(Z1);
R0=mvnrnd(mu0,sigma0,212);
R1=mvnrnd(mu1,sigma0*2,214);
figure;
subplot(2,2,1)
[p1]=anderson_2006_test(R0,R1,0.5,'on')

%%
%R0=MvnRnd_2(mu0,sigma0,212);
%R1=MvnRnd_2(mu1,sigma1,214);
R0=MvnRnd_2(mu0,sigma0,212);
R1=MvnRnd_2(mu1,sigma0*2,214);

subplot(2,2,2)
[p2]=anderson_2006_test(R0,R1,0.5,'on')


