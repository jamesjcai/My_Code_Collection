clear all; close all; clc;

  % space
  L=30; n=512;
  x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
  k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
  % time
  slices=40;
  t=linspace(0,2*pi,slices+1); dt=t(2)-t(1); 

  nf=1;
  slicesf=slices*nf;
  tf=linspace(0,2*pi*nf,slicesf+1);
  
 % initial conditions
  N=2;
  u=N*(sech(x)).';
  ut=fft(u);
  [t,utsol]=ode45('dmd_soliton_rhs',t,ut,[],k);
  for j=1:length(t)
      usol(j,:)=ifft(utsol(j,:));  % bring back to space
  end

subplot(2,2,1), waterfall(x,t,abs(usol)), colormap([0 0 0])
set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',[14],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('     x','Fontsize',[14]), ylabel('time','Fontsize',[14])
text(-40,8,3,'|u|','Fontsize',[14])
text(-20,3,7,'PDE','Fontsize',[14])

v = usol.';  % here is the data

%%%%%% body of DMD %%%%%%%%%%
v1 = v(:,1:end-1); v2 = v(:,2:end);

[U2,Sigma2,V2] = svd(v1, 'econ');
r=10; U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
Atilde = U'*v2*V/Sigma;
[W,D] = eig(Atilde);
Phi=v2*V/Sigma*W;

lambda=diag(D);
omega=log(lambda)/dt;

y0 = Phi\u;  % pseudo-inverse initial conditions

u_modes = zeros(r,length(t));  % DMD reconstruction for every time point
for iter = 1:length(tf)
    u_modes(:,iter) =(y0.*exp(omega*(tf(iter))));
end
u_dmd = Phi*u_modes;   % DMD resconstruction with all modes

figure(1), subplot(2,2,3), waterfall(x,1:r,abs(Phi).'), colormap([0 0 0])
set(gca,'Xlim',[-20 20],'Xtick',[-20 0 20],'Ylim',[1 r],'Ytick',[1 r],'Fontsize',[14],'Zlim',[0 0.6],'Ztick',[0 0.3 0.6])
xlabel('     x','Fontsize',[14]), ylabel('modes','Fontsize',[14])

subplot(2,2,4), plot(diag(Sigma2),'ko')


figure(1), subplot(2,2,2), waterfall(x,tf,abs(u_dmd).'), colormap([0 0 0])
   set(gca,'Ylim',[0 2*pi*nf],'Ytick',[0 3 6],'Fontsize',[14],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('     x','Fontsize',[14]), ylabel('time','Fontsize',[14])
text(-40,8,3,'|u|','Fontsize',[14])
text(-20,3,7,'DMD','Fontsize',[14])

for j=1:length(tf)
 error1(j)=norm(u_dmd(:,j)-usol(j,:).');
end


%% Extende DMD 
  clear all; close all; clc

  % space
  L=30; n=256;
  x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
  k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
  % time
  slices=2000;
  t=linspace(0,2*pi,slices+1); dt=t(2)-t(1); 

  nf=1;
  slicesf=slices*nf;
  tf=linspace(0,2*pi*nf,slicesf+1);
  
 % initial conditions
  N=2;
  u=N*(sech(x)).';
  ut=fft(u);
  [t,utsol]=ode45('dmd_soliton_rhs',t,ut,[],k);
  for j=1:length(t)
      usol(j,:)=ifft(utsol(j,:));  % bring back to space
  end

subplot(2,2,1), surfl(x,t,abs(usol)), shading interp, colormap(gray)
set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',[14],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('     x','Fontsize',[14]), ylabel('time','Fontsize',[14])
text(-40,8,3,'|u|','Fontsize',[14])
text(-20,3,7,'PDE','Fontsize',[14])

%% extended DMD
X=usol.'; X1=(X(:,1:end-1)); X2=(X(:,2:end));

tic
Ay1=X2*pinv(X1);
toc

tic
A1=X2*X1.';
A2=X1*X1.';
Ay2=A1*pinv(A2);
toc

figure(3), [w1,d1]=eig(Ay1); [w2,d2]=eig(Ay2); 
subplot(2,1,1), plot(log(diag(d1))/dt,'ro','Linewidth',[2]), hold on, plot(log(diag(d2))/dt,'ko','Linewidth',[2]), grid on
subplot(2,1,2), plot(log(diag(d1))/dt,'ro','Linewidth',[2]), hold on, plot(log(diag(d2))/dt,'ko','Linewidth',[2])
axis([-1000 500 -1000 1000])

%% kernel DMD 
  clear all; %close all; clc

  % space
  L=30; n=256;
  x2=linspace(-L/2,L/2,n+1); x=x2(1:n);
  k=(2*pi/L)*[0:n/2-1 -n/2:-1].';
  % time
  slices=200;
  t=linspace(0,2*pi,slices+1); dt=t(2)-t(1); 

  nf=1;
  slicesf=slices*nf;
  tf=linspace(0,2*pi*nf,slicesf+1);
  
 % initial conditions
  N=2;
  u=N*(sech(x)).';
  ut=fft(u);
  [t,utsol]=ode45('dmd_soliton_rhs',t,ut,[],k);
  for j=1:length(t)
      usol(j,:)=ifft(utsol(j,:));  % bring back to space
  end

figure(4), subplot(2,2,1), surfl(x,t,abs(usol)), shading interp, colormap(gray)
set(gca,'Ylim',[0 2*pi],'Ytick',[0 3 6],'Fontsize',[14],'Xlim',[-15 15],'Xtick',[-15 0 15],'Zlim',[0 4],'Ztick',[0 2 4])
xlabel('     x','Fontsize',[14]), ylabel('time','Fontsize',[14])
text(-40,8,3,'|u|','Fontsize',[14])
text(-20,3,7,'PDE','Fontsize',[14])


%% kernel
X=usol.'; X1=(X(:,1:end-1)); X2=(X(:,2:end));
[n,m]=size(X);

for j=1:m-1
 for jj=1:m-1
  YtYp(j,jj)=exp(-abs((X1(:,jj)-X2(:,j)).'*(X1(:,jj)-X2(:,j)))^2 );
   YsY(j,jj)=exp(-abs((X1(:,jj)-X2(:,j)).'*(X1(:,jj)-X2(:,j)))^2 );
 end 
end

[V,Sig]=eig(YsY);
Ay=(Sig*V')*(YtYp)*(V/Sig);

[W,D]=eig(Ay);
%plot(diag(D),'ko','Linewidth',[2])


figure(3), [w3,d3]=eig(Ay); 
subplot(2,1,1), plot(log(diag(d3))/dt,'bo','Linewidth',[2])
subplot(2,1,2), plot(log(diag(d3))/dt,'bo','Linewidth',[2]), axis([-1000 500 -1000 1000])

%% More kernels
p=20;
for j=1:m-1
 for jj=1:m-1       
  YtYp1(j,jj)=(1+((X1(:,jj)).')*(X2(:,j)))^p;
   YsY1(j,jj)=(1+((X1(:,jj)).')*(X1(:,j)))^p;
      
  YtYp2(j,jj)=(1+(abs(X1(:,jj)).')*abs(X2(:,j)))^p;
   YsY2(j,jj)=(1+(abs(X1(:,jj)).')*abs(X1(:,j)))^p;
      
  YtYp3(j,jj)=exp(-abs((X1(:,jj)-X2(:,j)).'*(X1(:,jj)-X2(:,j)))^2 );
   YsY3(j,jj)=exp(-abs((X1(:,jj)-X2(:,j)).'*(X1(:,jj)-X2(:,j)))^2 );
 
  YtYp4(j,jj)=exp(-abs((X1(:,jj).')*(X2(:,j))));
   YsY4(j,jj)=exp(-abs((X1(:,jj).')*(X2(:,j))));
 end 
end

figure(5)

[V,Sig]=eig(YsY1);
Ay=(Sig*V')*(YtYp1)*(V/Sig);
[W,D1]=eig(Ay);
subplot(2,2,1), plot(log(diag(D1))/dt,'ko','Linewidth',[2])
axis([-3000 3000 -100 100]), grid on, set(gca,'Fontsize',[14]), grid on

[V,Sig]=eig(YsY2);
Ay=(Sig*V')*(YtYp2)*(V/Sig);
[W,D2]=eig(Ay);
subplot(2,2,2), plot(log(diag(D2))/dt,'ko','Linewidth',[2])
axis([-3000 3000 -100 100]), grid on,set(gca,'Fontsize',[14]), grid on

[V,Sig]=eig(YsY3);
Ay=(Sig*V')*(YtYp3)*(V/Sig);
[W,D3]=eig(Ay);
subplot(2,2,3), plot(log(diag(D3))/dt,'ko','Linewidth',[2])
axis([-3000 3000 -100 100]), grid on, set(gca,'Fontsize',[14]), grid on

[V,Sig]=eig(YsY4);
Ay=(Sig*V')*(YtYp4)*(V/Sig);
[W,D4]=eig(Ay);
subplot(2,2,4), plot(log(diag(D4))/dt,'ko','Linewidth',[2])
axis([-3000 3000 -100 100]), grid on, set(gca,'Fontsize',[14]), grid on



figure(3),
subplot(2,1,1), set(gca,'Fontsize',[14]), grid on
subplot(2,1,2), set(gca,'Fontsize',[14]), grid on
%legend('     ','     ','     ','Fontsize',[15])





  
  