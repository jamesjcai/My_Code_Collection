n=1024;
[I,J]=meshgrid(1:n,1:n);
DFT=w.^((I-1).*(J-1));
imagesc(real(DFT))
issymmetric(DFT)


%%

dt=0.001;
t=0:dt:1;
fclean=sin(2*pi*50*t)+sin(2*pi*120*t);
f=fclean+2.5*randn(size(t));
figure; 
subplot(3,1,1)
plot(t,f); hold on;
plot(t,fclean);
subplot(3,1,3)
plot(1:length(f),abs(fft(f)/length(f)));

