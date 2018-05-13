%rng('default');
load('testdata_scz_Z1_ctl_Z0.mat')
mu0=mean(Z0);
sigma0=cov(Z0);
mu1=mean(Z1);
sigma1=cov(Z1);

%%
smpsiz=100:100:2100;
lamdav=1:0.2:5;
Zc=cell(21,21);
%%
for kx=1:21
    for ky=1:21
        Sz=smpsiz(kx);
        Lv=linspace(1,lamdav(ky),round(Sz/4));
        p=ones(100,1);
        [kx ky]
        parfor k=1:100
            R0=mvnrnd(mu0,sigma0,round(Sz/2));
            R1a=mvnrnd(mu1,sigma1,round(Sz/4));
            R1b=[];
            for k2=1:round(Sz/4)
                R1b=[R1b;mvnrnd(mu1,sigma1*Lv(k2),1)];
            end
            R1=[R1a;R1b];
            %subplot(2,2,2)
            p(k)=anderson_2006_test(R0,R1,0.5,'off');
        end
        Zc{kx,ky}=p;
    end
end
[X,Y] = meshgrid(smpsiz,lamdav);
save zxy2 Zc X Y
%%

Z=cellfun(@(X)(sum(X<0.01)),Zc);
figure; mesh(X,Y,Z)
figure; contour(X,Y,Z,'showtext','on')

