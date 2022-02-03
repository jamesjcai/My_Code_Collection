M = 15;  % number of measurements

% projType = 1; % uniform random projection
projType = 2; % Gaussian random projection
% projType = 3; % Single pixel measurement

%% Create random measurement matrix
C = zeros(M,n*n);
Theta = zeros(M,n*n);
xmeas= zeros(n,n);
for i=1:M
    xmeas = 0*xmeas;
    if(projType==1)
        xmeas = rand(n,n);
    elseif(projType==2)
        xmeas = randn(n,n);
    elseif(projType==3)
        xmeas(ceil(n*rand),ceil(n*rand)) = 1;
    end
    C(i,:) = reshape(xmeas,n*n,1);
    Theta(i,:) = reshape((ifft2(xmeas)),1,n*n);
end

%% Project data
YDAT = C*XDAT;
if(saveFLAG)
    save([filename,'_DATAY.mat']);
end

%% plot C
figure, colormap bone
Cimg = sum(C,1);
Ximg = XDAT(:,1).*Cimg';
imagesc(reshape(Ximg,n,n));
figure, colormap bone;
surf(X,Y,Z,reshape(Cimg,n,n));
view(10,32)
axis equal, axis off
drawnow