%rng('default');
load('testdata_scz_Z1_ctl_Z0.mat')
mu0=mean(Z0);
sigma0=cov(Z0);
mu1=mean(Z1);
sigma1=cov(Z1);

%%
smpsiz=400;
Sz=smpsiz;
Lamda=3.5;
lamdav=linspace(1,Lamda,round(Sz/4));
%%
R0=mvnrnd(mu0,sigma0,round(Sz/2));
R1a=mvnrnd(mu1,sigma1,round(Sz/4));
R1b=[];
for k=1:round(Sz/4)
    R1b=[R1b;mvnrnd(mu1,sigma1*lamdav(k),1)];
end
R1=[R1a;R1b];
figure;
subplot(2,2,1)
p1=anderson_2006_test(R0,R1,0.5,'on')
ylim([0 19])
subplot(2,2,2)
p2=anderson_2006_test(Z0,Z1,0.5,'on')
