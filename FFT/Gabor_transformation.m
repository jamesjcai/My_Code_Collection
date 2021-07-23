% https://www.mathworks.com/matlabcentral/answers/757959-gabor-transform-for-a-sine-wave
% https://youtu.be/EfWnEldTyPA

fy=5; %signal frequency in Hz
wy=2*pi*fy; %signal frequency in rad/s
fs=100; %sampling frequency in Hz
tiv=1/fs; %time interval between samples;
t=0:tiv:10; %time intervals 
y=sin(wy*t); %signal data set
subplot(311),plot(t,y); %plots figure
axis([0 10 -1.5 1.5]);
xlabel('seconds'); title('sine signal');
subplot(312),stem(t,y,'*'); %plots figure
axis([0 1 -1.5 1.5]);
xlabel('seconds'); title('Samples of sine signal for Time Period of 1 second');
subplot(313),cqt(y,'SamplingFrequency',fs);%plots histogram of Constant-Q nonstationary Gabor transform

%%

fy=5; %signal frequency in Hz
wy=2*pi*fy; %signal frequency in rad/s
fs=100; %sampling frequency in Hz
tiv=1/fs; %time interval between samples;
t=0:tiv:10; %time intervals 
y=sin(wy*t); %signal data set
subplot(311),plot(t,y); %plots figure
axis([0 10 -1.5 1.5]);
xlabel('seconds'); title('sine signal');
subplot(312),stem(t,y,'*'); %plots figure
axis([0 1 -1.5 1.5]);
xlabel('seconds'); title('Samples of sine signal for Time Period of 1 second');
subplot(313),cqt(y,'SamplingFrequency',fs);%plots his

