%% make the example movie
% parameters
nx = 80; % horizontal pixels
ny = 80; % vertical pixels
n  = nx * ny;
T = 10; % sec
dt = 0.01;
t = dt:dt:T;
sig = 0.1; % magnitude of gaussian noise added

% 3 underlying modes
[Xgrid, Ygrid] = meshgrid(1:nx, 1:ny);

mode1.u = exp(-((Xgrid-40).^2/250+(Ygrid-40).^2/250));
mode1.f = 5.55; % cycles/sec
mode1.A = 1;
mode1.lambda = exp(1i * mode1.f*2*pi*dt);
mode1.range = [0, 5];

mode2.u = zeros(nx, ny); 
mode2.u(nx-40:nx-10, ny-40:ny-10) = 1;
mode2.f = 0.9; % cycles/sec 
mode2.A = 1;
mode2.lambda = exp(1i * mode2.f*2*pi*dt);
mode2.range = [3, 7];

mode3.u = zeros(nx, ny); 
mode3.u(1:nx-20, 1:ny-20) = 1;
mode3.f = 0.15; % cycles/sec 
mode3.A = 0.5;
mode3.lambda = exp(1i * mode3.f*2*pi*dt);
mode3.range = [0, T];

modes(1) = mode1;
modes(2) = mode2;
modes(3) = mode3;

% make the movie
Xclean = zeros(n, numel(T));
for ti = 1:numel(t),
    Snap = zeros(nx, ny);
    for mi = 1:numel(modes),
        mymode = modes(mi);
        if ti > round(mymode.range(1)/dt) && ...
            ti < round(mymode.range(2)/dt),
            Snap = Snap + mymode.A * mymode.u * ...
                real(mymode.lambda^ti);
        end;
    end;
    Xclean(:, ti) = Snap(:);
end;

% add noise
Noise = sig * randn(size(Xclean));
X = Xclean + Noise;

figure; 
imagesc(X);

%% take a look at the movie
figure;
for ti = 1:numel(t),
    Pic = reshape(X(:, ti), nx, ny);
    imagesc(Pic, range(X(:))/2*[-1 1]);
    axis square; axis off;
    pause(dt);
end;

%% compute mrDMD
L = 6; % number of levels
r = 10; % rank of truncation

mrdmd = mrDMD(X, dt, r, 2, L);

% compile visualization of multi-res mode amplitudes
[map, low_f] = mrDMD_map(mrdmd);
[L, J] = size(mrdmd);

%%
figure; 
imagesc(-map); 
set(gca, 'YTick', 0.5:(L+0.5), 'YTickLabel', floor(low_f*10)/10); 
set(gca, 'XTick', J/T*(0:T) + 0.5);
set(gca, 'XTickLabel', (get(gca, 'XTick')-0.5)/J*T);
axis xy;
xlabel('Time (sec)');
ylabel('Freq. (Hz)');
colormap pink;
grid on;

%%
figure;
imagesc(reshape(abs(mrdmd{1,1}.Phi(1:n,1)), nx, ny));
axis square;