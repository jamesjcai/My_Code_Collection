clear; clc; close all;

%
%This demonstrates the usage of the Dynamic Mode Decomposition wrapper for
%a shedding cylinder. Work done for Dr. Louis Cattafesta's Flow Control
%class at Florida State University - data set and directions provided by him. 
%Wrapper built to facilitate data processing of other flow fields.


%%
%==============PERFORMS THE DMD AND DATA RESHAPING==============
load('VelocityFieldData.mat' );
dt=0.3638; %Delta t given [s]

Vall(:,:,:,1)=permute(vx,[3 1 2]);
Vall(:,:,:,2)=permute(vy,[3 1 2]);
r=25;

[Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, GrowthRates, POD_Mode_Energies]=dmd_rom(Vall, r, dt);


%Finds dominant oscillatory mode
dominantModeNo=find(ModeAmplitudes==max(ModeAmplitudes));

%%
%==========Mode Amplitudes===========
figure('color', 'w');
stem(abs(ModeFrequencies),abs(ModeAmplitudes), 'k^', 'filled'); hold on
plot(abs(ModeFrequencies(dominantModeNo)),abs(ModeAmplitudes(dominantModeNo)),'ro','MarkerSize',10);
set(gca, 'YScale', 'log')
title('Mode amplitude versus Frequency');
xlabel('f [Hz]');ylabel('Amplitude');



%%
%==========Plot Mode Shape===========
m=1;
skip=4;

theta=linspace(0,6*pi,100);

figure('color', 'w'); 
for i=1:length(theta)
    U=squeeze(real(Eigenvectors(m,:,:,1)*exp(1i*theta(i))));
    V=squeeze(real(Eigenvectors(m,:,:,2)*exp(1i*theta(i))));
    
    q=quiver(U(1:skip:end,1:skip:end),V(1:skip:end,1:skip:end),4);
    set(q,'Color','k')
    drawnow;
    pause(0.03);
end
