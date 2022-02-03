%% Define time and space discretizations
n = 200;
m = 80;
x = linspace(-15,15,n);
t = linspace(0,8*pi,m); 
dt = t(2) - t(1); 
[Xgrid,T] = meshgrid(x,t);

%% Create two spatio-temporal patterns
f1 = 0.5*cos(Xgrid) .* (1+0*T);  % time-independent!
f2 = (sech(Xgrid).*tanh(Xgrid)) .* (2*exp(1j*2.8*T));

%% Combine signals and make data matrix
X = (f1 + f2)'; % Data Matrix

figure;
surfl(real(X')); 
shading interp; colormap(gray); view(-20,60);