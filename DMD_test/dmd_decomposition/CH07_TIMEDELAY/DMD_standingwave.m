clear all, close all, clc

dt = 0.01;
t = 0:dt:10;
x = sin(t);
plot(t,x), hold on

X = x(1:end-1);
X2 = x(2:end);

[U,S,V] = svd(X,'econ');
Atilde = U'*X2*V*inv(S);
[W,Lambda] = eig(Atilde);
Omega = log(Lambda)/dt
phi = X2*V*inv(S)*W;

b = (phi*exp(Omega*t))'\x';
xdmd = b*phi*exp(Omega*t);
plot(t,xdmd,'r--')

%% Augmented DMD
Xaug = [x(1:end-2);
    x(2:end-1)];
Xaug2 = [x(2:end-1);
    x(3:end)];

[U,S,V] = svd(Xaug,'econ');
Atilde = U'*Xaug2*V*inv(S);
[W,Lambda] = eig(Atilde);
Omega = diag(log(diag(Lambda)))/dt
Phi = Xaug2*V*inv(S)*W;

b = Phi\Xaug(:,1);
for k=1:length(t)
    xaugdmd(:,k) = Phi*exp(Omega*t(k))*b;
end
plot(t,real(xaugdmd(1,:)),'b--','LineWidth',1.5)
ylim([-1 1]), grid on
legend('Data','DMD','Augmented DMD')
xlabel('Time'), ylabel('State')


set(gcf,'Position',[100 100 500 200])
set(gcf,'PaperPositionMode','auto')
% print('-depsc2', '-loose', 'standingwave');