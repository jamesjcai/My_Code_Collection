% cook up linear dynamical system
% compare DMD with forward/backward DMD and total least square DMD
clear all, close all, clc    

%%
figure_dir = '../figures/';
show_legend = 0;

% parameters
sig = 0.5;  % noise standard deviation
niterations = 500;  % size of random ensemble
m = 100;  % number of snapshots


A = [1, 1; -1, 2]/sqrt(3);
n = 2;

X = zeros(n, m);
X(:, 1) = [0.5, 1]';

for k = 2:m,
    X(:, k) = A * X(:, k-1);
end;

lambda = eig(A);

%%
mylambdas = [];
for iter = 1:niterations,

    Y = X + sig*randn(size(X));
    
    % build data matrices
    Y1 = Y(:, 1:end-1);
    Y2 = Y(:, 2:end);

    % compute DMD
    [U, S, V] = svd(Y1, 'econ');
    Atilde = U' * Y2 * V / S;

    [W, D] = eig(Atilde);
    lam = diag(D); % eigenvalues
    
    mylambdas = [mylambdas; lam];
end;


%% now forward/backward dmd

fblambdas = [];
for iter = 1:niterations,

    Y = X + sig*randn(size(X));
    
    % build data matrices
    Y1 = Y(:, 1:end-1);
    Y2 = Y(:, 2:end);

    % compute forward DMD
    [U, S, V] = svd(Y1, 'econ');
    f_Atilde = U' * Y2 * V / S;
    
    % compute backward DMD
    [U, S, V] = svd(Y2, 'econ');
    b_Atilde = U' * Y1 * V / S;
    
    % estimate forward/backward DMD
    Atilde = (f_Atilde * inv(b_Atilde)) ^ 0.5;
    [W, D] = eig(Atilde);
    lam = diag(D);
    
    fblambdas = [fblambdas; lam];
end;


%%
neigh = sig * [-1 0.5];
gray = 0.7 * [1 1 1];

figure;
hold on;
plot(real(mylambdas), imag(mylambdas), '.', 'Color', 'b'); 
plot(real(fblambdas), imag(fblambdas), '+', 'Color', gray); 
plot(real(lambda), imag(lambda), 'r^', 'MarkerFaceColor', 'r');
rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
    'LineStyle', '--');
axis square;
l1 = legend('DMD', 'fbDMD', 'true');
set(l1,'FontSize',14)
xlim(real(lambda(1)) + neigh);
ylim(imag(lambda(1)) + neigh);
title(sprintf('Noise = %g', sig));

%% total least square DMD

r = 2;

tlslambdas = [];
for iter = 1:niterations,

    Y = X + sig*randn(size(X));
    
   

    % build data matrices
    Y1 = Y(:, 1:end-1);
    Y2 = Y(:, 2:end);

%     [U_Y, S_Y, V_Y] = svd(Y1, 'econ');
%     Y1tilde = U_Y(:, 1:r)' * Y1;
%     
%     [U_Y, S_Y, V_Y] = svd(Y2, 'econ');
%     Y2tilde = U_Y(:, 1:r)' * Y2;
%     
%     Z = [Y1tilde; Y2tilde];
     
    % and stack
    Z = [Y1; Y2];
    
    % tls DMD
    [U, S, V] = svd(Z, 'econ');
    U11 = U(1:r, 1:r);
    U21 = U(r+1:end, 1:r);
    
    Atilde = U21 * inv(U11);
    
    [W, D] = eig(Atilde);
    lam = diag(D);
    
    tlslambdas = [tlslambdas; lam];
end;


%%
neigh = sig * [-0.9 0.5];
gray = 0.7 * [1 1 1];

figure;
hold on;
rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
    'LineStyle', '--');
plot(real(mylambdas), imag(mylambdas), '.', 'Color', 'b'); 
plot(real(tlslambdas), imag(tlslambdas), 'o', 'Color', gray); 
plot(real(fblambdas), imag(fblambdas), '+', 'Color', 'k'); 
plot(real(lambda), imag(lambda), 'r^', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 8);
axis square;
if show_legend,
    legend({'DMD', 'tlsDMD', 'fbDMD', 'true'}, ...
        'Location', 'best');
end;
xlim([0.4 1.1]);
ylim([0 0.8]);
% set(gca,'XTick',[.4 .6 .8 1],'YTick',[0 .2 .4 .6 .8])
l1 = legend('DMD', 'tlsDMD','fbDMD', 'true');
% set(l1,'FontSize',14)
xlim(real(lambda(1)) + neigh);
ylim(imag(lambda(1)) + neigh);
% title(sprintf('m = %i, sigma = %g', m, sig));


set(gcf, 'Color', 'w', 'Position', [500 500 200 200]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4 4], 'PaperPositionMode', 'manual');
% print('-depsc', '-loose', [figure_dir sprintf('DMD_denoise_zoom_m%i_sigma%g.eps', m, sig)]);