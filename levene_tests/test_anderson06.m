rng('default');
load('testdata_scz_Z1_ctl_Z0.mat')
mu0=mean(Z0);
sigma0=cov(Z0);
mu1=mean(Z1);
sigma1=cov(Z1);

%%
lamdav=1:0.1:1.8;
figure;
for kx=1:9
    Lv=lamdav(kx);
    p=[];
    parfor k=1:100
        k
        R0=mvnrnd(mu0,sigma0,212);
        R1a=mvnrnd(mu1,sigma1,107);
        R1b=mvnrnd(mu1,sigma1,107);
        R1=[R1a;R1b];
        %figure;
        %subplot(2,2,1)
        [p1]=anderson_2006_test(R0,R1,0.5,'off');

        R0=mvnrnd(mu0,sigma0,212);
        R1a=mvnrnd(mu1,sigma1,107);
        R1b=mvnrnd(mu1,sigma1*Lv,107);
        R1=[R1a;R1b];
        %subplot(2,2,2)
        [p2]=anderson_2006_test(R0,R1,0.5,'off');
        p=[p;[p1 p2]];
    end
    p1=p(:,1);
    p2=p(:,2);    
    subplot(3,3,kx)
    qqplot(-log10(p1),-log10(p2))
    refline([1 0])
    box on
    title(sprintf('lamda=%.2f',Lv));
end
