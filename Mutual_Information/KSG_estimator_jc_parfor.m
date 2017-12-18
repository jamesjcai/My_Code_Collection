function [I1, I2] = KSG_estimator_jc_parfor(X, Y, k)
% tic

N = size(X, 1);

% compute distance between each sample and its k-th nearest neighbour
dz = zeros(N, N);
%dx = zeros(N, N);
%dy = zeros(N, N);

dx=squareform(pdist(X,'euclidean'));
dy=squareform(pdist(Y,'euclidean'));
%  for i = 1:N
%      for j = 1:N
% % %        dx(i,j) = sqrt(sum((X(i, :) - X(j, :)).^2));
% % %        dy(i,j) = sqrt(sum((Y(i, :) - Y(j, :)).^2));
%          dz(i,j) = max([dx(i, j), dy(i, j)]);
%      end
%  end
dff=dx-dy;
dz(dff>=0)=dx(dff>=0);
dz(dff<0)=dy(dff<0);



% dv = abs(bsxfun(@minus,v,v'));
% toc 
% 
% tic

% find nx(i) and ny(i)
Eps = zeros(N, 1);
Nn = zeros(N, 1);

nx1 = zeros(N, 1);
ny1 = zeros(N, 1);
nx2 = zeros(N, 1);
ny2 = zeros(N, 1);


for i = 1:N
    
    curr_dx_col = dx(i, :);
    curr_dx_col(i) = [];
    
    curr_dy_col = dy(i, :);
    curr_dy_col(i) = [];
    
    dz_col = dz(i, :);
    dz_col(i) = [];
    
    [EpsSample, NnSample] = sort(dz_col, 'ascend');
    
   
    Eps(i) = EpsSample(k);
    Nn(i) = NnSample(k);
    
    nx1(i) = sum(curr_dx_col < Eps(i));
    ny1(i) = sum(curr_dy_col < Eps(i));
    
    nx2(i) = sum(curr_dx_col <= Eps(i));
    ny2(i) = sum(curr_dy_col <= Eps(i));
    
end


nObs=N;
% mutual information estimators
I1 = psi(k) - sum(psi(nx1 + 1) + psi(ny1 + 1)) / N + psi(N);

if nargout>1
I2 = psi(k) - 1/k - sum(psi(nx2) + psi(ny2)) / nObs + psi(nObs);
end

% if I1 < 0
%     warning('First estimator is negative -> 0');
%     I1 = 0;
% end
% 
% if I2 < 0
%     warning('Second estimator is negative -> 0');
%     I2 = 0;
% end

% toc 
end
