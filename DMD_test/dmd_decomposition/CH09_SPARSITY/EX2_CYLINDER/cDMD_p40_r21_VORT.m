clear all, close all, clc
load ../../DATA/FLUIDS/CYLINDER_ALL.mat
p = 40;
X = VORTALL;
X2 = [VORTALL(:,2:end) VORTEXTRA];
[U,S,V] = svd(X,'econ');

% projType = 1; % uniform random projection
projType = 2; % Gaussian random projection
% projType = 3; % Single pixel measurement

C = zeros(p,m*n);
Theta = zeros(p,m*n);
xmeas= zeros(m,n);
for i=1:p
    i
    xmeas = 0*xmeas;
    if(projType==1)
        xmeas = rand(m,n);
    elseif(projType==2)
        xmeas = randn(m,n);
    elseif(projType==3)
        xmeas(ceil(m*rand),ceil(n*rand)) = 1;
    end
    C(i,:) = reshape(xmeas,n*m,1);
    Theta(i,:) = reshape((ifft2(xmeas)),1,n*m);
end

Y = C*X;
Y2 = C*X2;
[UY,SY,VY] = svd(Y,'econ');

%% figure, plot singular values
figure
semilogy(diag(S(1:20,1:20)),'-ok','LineWidth',1.5)
hold on
semilogy(diag(SY(1:20,1:20)),'or','LineWidth',1.5)
xlim([0 20])
ylim([2 4000])
grid on
set(gca,'YTick',[10 100 1000],'YMinorGrid','on','YMinorTick','on')
set(gcf,'Position',[100 100 150 130])
set(gcf,'PaperPositionMode','auto')
fname2 = 'figures/csDMD_p600_T2_r21/csDMD_VORTdiagS.eps';
% print('-depsc2', '-loose', fname2);


figure
plot(cumsum(diag(S(1:20,1:20)))/sum(diag(S)),'-ok','LineWidth',1.5)
hold on
plot(cumsum(diag(SY(1:20,1:20)))/sum(diag(SY)),'or','LineWidth',1.5)
xlim([0 20])
grid on
set(gcf,'Position',[100 100 150 130])
set(gcf,'PaperPositionMode','auto')
fname2 = 'figures/csDMD_p600_T2_r21/csDMD_VORTcumsumdiagS.eps';
% print('-depsc2', '-loose', fname2);

%%
r = 21;
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);

Phi = X2*V*inv(S)*W;

UY = UY(:,1:r);
SY = SY(1:r,1:r);
VY = VY(:,1:r);
AtildeY = UY'*Y2*VY*inv(SY);
[WY,eigsY] = eig(AtildeY);

PhiY = Y2*VY*inv(SY)*WY;
PhiXtilde = X2*VY*inv(SY)*WY;

%% COMPRESSED MODES FROM PHIXTILDE
for i=1:r
    plotCylinderNoSave(reshape(real(PhiXtilde(:,i)),nx,ny));
    plotCylinderNoSave(reshape(imag(PhiXtilde(:,i)),nx,ny));
end