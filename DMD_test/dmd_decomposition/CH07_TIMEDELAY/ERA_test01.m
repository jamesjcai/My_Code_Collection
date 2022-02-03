clear all, close all, clc

dt = 0.01;
t = 0:dt:10;
x = exp(-t) + exp(-3*t);
plot(t,x)

H = [x(1:5);
    x(2:6);
    x(3:7);
    x(4:8);
    x(5:9)];

[U,S,V] = svd(H);
r = rank(S)

%%
alpha = -.01;
beta = -10;
t = 0:.01:100;
x = exp(-.01*t)+exp(-10*t);
plot(t,x,'k'), hold on

H = [x(1:end-2);
    x(2:end-1)];
H2 = [x(2:end-1);
    x(3:end)];
[U,S,V] = svd(H,'econ');
r = rank(S)
YY(1,1,:) = x;
[Ar,Br,Cr,Dr,HSVs] = ERA(YY,100,100,1,1,r);
sys = ss(Ar,Br,Cr,Dr,dt);
u = zeros(size(t));
u(1) = 1;
[y,t] = lsim(sys,u,t); % discrete impulse
plot(t,y,'--r')
xlabel('Time'), ylabel('State')
legend('Data','ERA Model')
log(eig(Ar))/dt

%%
figure
subplot(1,2,1)
plot(t,x,'k'), hold on
plot(t,y,'r--');
xlabel('Time')
ylabel('State')
xlim([0 1])
legend('Data','ERA Model')
subplot(1,2,2)
plot(t,x,'k'), hold on
plot(t,y,'r--');
xlabel('Time')
xlim([0 100])
set(gcf,'Position',[100 100 500 175])


set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', 'ERAtest01');