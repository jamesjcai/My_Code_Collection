clc; clear ;close all;
 x = pi/100:pi/100:2*pi;T=200;
  gam1 = rgam(T, 3, 1, 1);
   y = sin(x); 
    X1=[x;y];
 x=interp1(linspace(0,1,200),x,gam1,'spline');
 y=interp1(linspace(0,1,200),y,gam1,'spline');
 X2=[x;y];
% % % X1=y;
% % T=200;

%  T=length(X1);
%  load toydata.mat;X1=C(:,:,30);
%   T=100;
% gam1 = rgam(T, 3, 1, 1);
% for  i=1:T-1
%     tmp=gam1(i)*(T-1)+1;
%     t0=floor(tmp);
% %     if t0<1
% %         tmp=1;
% %     else
%     t1=ceil(tmp);
%     p0=X1(t0);
%     p1=X1(t1);
%     if t0==t1;
%        X2(i)=p0;
%     else                
%        X2(i)=(p1-p0)/(t1-t0)*(tmp-t0)+p0;
%     end
% %     end
% end
%    X2(T)=X1(T);
%  X2=interp1(linspace(0,1,T),X1,gam1);
%   X1=curve_to_q(X1);
%  X2 = Group_Action_by_Gamma_Coord(X1,gam1);
%     X2=q_to_curve(X2);
save X1;save X2;save gam1;
clear;

%%
% load X1;
load X2;

[dist,X2n,q2n,X1x,q1,gam]=mygeod(X1,X2);

%%
figure;
plot(linspace(0,1,T),gam1,'k');hold on;plot(linspace(0,1,200),gam,'r');

