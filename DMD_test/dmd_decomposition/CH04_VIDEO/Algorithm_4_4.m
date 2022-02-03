%% Compute DMD Foreground Solution
b = Phi_bg \ X(:, 1);
X_bg = zeros(numel(omega_bg), length(t));
for tt = 1:length(t),
    X_bg(:, tt) = b .* exp(omega_bg .* t(tt));
end;
X_bg = Phi_bg * X_bg;
X_bg = X_bg(1:n, :);

figure;
surfl(real(X_bg')); 
shading interp; colormap(gray); view(-20,60);

%% Compute DMD Background Solution
b = Phi_fg \ X(:, 1);
X_fg = zeros(numel(omega_fg), length(t));
for tt = 1:length(t),
    X_fg(:, tt) = b .* exp(omega_fg .* t(tt));
end;
X_fg = Phi_fg * X_fg;
X_fg = X_fg(1:n, :);

figure;
surfl(real(X_fg')); 
shading interp; colormap(gray); view(-20,60);