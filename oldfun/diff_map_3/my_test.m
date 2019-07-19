load MY_data.mat

data=X';
thr=15;
genes=genelist;
L=25;

censoredornot=0; %=1 if you want censoring otherwise set it to 0
t=1;
msteps=10;
stepsize=0.2;

% %%%%%%%%%%%%%%%this part is to estimate a proper sigma, not needed if a suitable sigma is already deceided%%%
% %with censoring this part can take about 7 mins on typical PCs. Without
% %censoring much faster (about one second)
% begin=0;
% [logsigma,dim_norm] = diffusion_map_kt(data,censoredornot,thr,thr,thr+L,begin,msteps,stepsize);
% [dim_max,id]=max(dim_norm); %find local maxima

[psi,E] = diffusion_map_main(data,censoredornot,thr,thr,thr+L, t, sigma);
%%%