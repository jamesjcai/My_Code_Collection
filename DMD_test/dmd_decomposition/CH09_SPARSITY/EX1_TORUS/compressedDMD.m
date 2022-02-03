%% COMPUTE DMD
% step 0
YD1 = YDAT(:,1:end-1);
YD2 = YDAT(:,2:end);
% step 1
[U,S,V] = svd(YD1,0);
% step 2
r=10; % for two mode and K=5
Sinv = S(1:r,1:r)^(-1);
Atilde = U(:,1:r)'*YD2*V(:,1:r)*Sinv(1:r,1:r);
% step 3
[W,D] = eig(Atilde);
% step 4
PhiY = YD2*V(:,1:r)*Sinv*W;
PhiXr = XD2*V(:,1:r)*Sinv*W;  % reconstructed

%% Plot DMD eigenvalues
figure, hold on
dmdeigsFULL = dmdeigs;
dmdeigs = log(eig(Atilde))*100;
lambdap = damping + sqrt(-1)*F*2*pi;
lambdan = damping - sqrt(-1)*F*2*pi;
lambda = [lambdap lambdan];
scatter(real(lambda),imag(lambda),'kx')
scatter(real(dmdeigsFULL),imag(dmdeigsFULL),'bo')
scatter(real(dmdeigs),imag(dmdeigs),'rd')
box on, grid on
legend('True eigenvalues','Exact DMD eigenvalues','Compressed DMD eigenvalues')

%% Plot Reconstructed DMD modes from X
figure
for k=1:K
    for j=1:2
        subplot(2,K,(j-1)*K+k)
        if(j==1)
            ximg = real(reshape(PhiXr(:,(k-1)*2+1)+PhiXr(:,k*2),n,n));
        end
        if(j==2)
            ximg = imag(reshape(PhiXr(:,(k-1)*2+1)-PhiXr(:,k*2),n,n));
        end
        surf(X,Y,Z,(ximg),'EdgeColor','none');
        view(10,32), colormap jet, axis equal, axis off, drawnow
    end
end

%% Reconstruct Modes with Compressive Sensing.
for i=1:K
    y = PhiY(:,2*i-1);    
%     cvx_begin;
%         variable ahat(n*n,1);
%         minimize( norm (ahat,1) );
%         subject to
%         Theta*ahat == y;
%     cvx_end;
    ahat = cosamp(Theta,y,2,10^-5,100);
    Phi(:,2*i-1) = reshape(real(ifft2(reshape(ahat,n,n))),n*n,1);
    Phi(:,2*i) = reshape(imag(ifft2(reshape(ahat,n,n))),n*n,1);
end
%% Plot Modes
figure
for k=1:K
    for j=1:2
        subplot(2,K,(j-1)*K+k)
        if(j==1)
            ximg = reshape(Phi(:,(k-1)*2+1),n,n);
        end
        if(j==2)
            ximg = reshape(Phi(:,k*2),n,n);
        end
        MODESCSDMD(j,k,:) = ximg(:);

        surf(X,Y,Z,(ximg),'EdgeColor','none');
        view(10,32), colormap jet, axis equal, axis off
    end
end