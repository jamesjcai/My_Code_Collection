addpath ./utils/
%% Generate signal
n = 4096, p = 256;
t = linspace(0, 1, n);
x = sin(73*2*pi*t) + sin(531*2*pi*t);

%% Randomly sample signal
perm = round(rand(p, 1) * n);
y = x(perm)';

%% Form matrix operators
Psi = dct(eye(n, n));
CPsi = Psi(perm, :);

%% L1 minimization (through linear program)
s = cosamp(CPsi,y,10,1.e-10,10);
xreconstruct = idct(s);

%% 
plot(t,x)
hold on
plot(t,xreconstruct)
xlim([.4 .5])