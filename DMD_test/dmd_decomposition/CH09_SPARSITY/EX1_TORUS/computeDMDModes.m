%% COMPUTE DMD
% step 0
XD1 = XDAT(:,1:end-1);
XD2 = XDAT(:,2:end);
% step 1
[U,S,V] = svd(XD1,0);
% step 2
Sinv = S(1:r,1:r)^(-1);
Atilde = U(:,1:r)'*XD2*V(:,1:r)*Sinv(1:r,1:r);
% step 3
[W,D] = eig(Atilde);
% step 4
Phi = XD2*V(:,1:r)*Sinv*W;

%% Plot DMD eigenvalues
figure, hold on
dmdeigs = log(eig(Atilde))*100;
lambdap = damping + sqrt(-1)*F*2*pi;
lambdan = damping - sqrt(-1)*F*2*pi;
lambda = [lambdap lambdan];
scatter(real(lambda),imag(lambda),'kx')
scatter(real(dmdeigs),imag(dmdeigs),'rd')  
box on, grid on
legend('True eigenvalues','DMD eigenvalues')
set(gcf,'Position',[100 100 400 300])

%% Plot Modes
figure
for k=1:K
    for j=1:2
        subplot(2,K,(j-1)*K+k)
        if(j==1)
            ximg = real(reshape(Phi(:,(k-1)*2+1)+Phi(:,k*2),n,n));
        end
        if(j==2)
            ximg = imag(reshape(Phi(:,(k-1)*2+1)-Phi(:,k*2),n,n));
        end
        MODESDMD(j,k,:) = ximg(:);
        surf(X,Y,Z,(ximg),'EdgeColor','none');        
        view(10,32), colormap jet, axis equal, axis off
    end
end