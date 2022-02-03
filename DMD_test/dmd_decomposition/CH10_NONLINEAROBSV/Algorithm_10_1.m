%% Define spatial discretization
L=30; n=512;  %% Domain length and # points
xi2=linspace(-L/2,L/2,n+1);  %% domain discretization
xi=xi2(1:n);  %% periodic domain
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';  %%  wavenumbers

%% Define time discretization
slices=20; 
t=linspace(0,pi,slices+1); dt=t(2)-t(1);

%% Create initial conditions
q=2*(sech(xi)).';
qt=fft(q);
%% Combine signals

%% Solve with Runge-Kutta
[t,qtsol]=ode45('dmd_soliton_rhs',t,qt,[],k);

%% Bring back to time domain and store data
for j=1:length(t)
   qsol(j,:)=ifft(qtsol(j,:));
end
X = qsol.';