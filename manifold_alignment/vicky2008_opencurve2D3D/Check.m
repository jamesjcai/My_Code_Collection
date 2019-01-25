clc; clear ;close all;
x = 0:pi/100:2*pi;
y = sin(x);
X1=[x;y];
T=201;
[q1] = curve_to_q(X1);
gam1 = rgam(T, 3, 1, 1);
X2 = Group_Action_by_Gamma_Coord(q1,gam1);
X2 = q_to_curve(X2);
save X1;save X2;save gam1;
clear;
load X1;load X2;
lam=0;
[q1] = curve_to_q(X1);
[q2] = curve_to_q(X2);
[gam] = DynamicProgrammingQ(q2/sqrt(InnerProd_Q(q2,q2)),q1/sqrt(InnerProd_Q(q1,q1)),lam,0);
 gamI = invertGamma(gam);
 gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    X2n = Group_Action_by_Gamma_Coord(X2,gam);
    q2n = curve_to_q(X2n);
    figure(1);
    plot(linspace(0,1,T),gam1,'k');hold on;plot(linspace(0,1,T),gamI,'r');
%  
% % Resample the curves to have N points
%     N = 200;
%     X1 = ReSampleCurve(X1,N);
%     X2 = ReSampleCurve(X2,N);
% 
% %Center curves, not really needed but good for display purposes
%     X1 = X1 - repmat(mean(X1')',1,size(X1,2));
%     X2 = X2 - repmat(mean(X2')',1,size(X2,2));    
%     
% % Form the q function for representing curves and find best rotation
%     [q1] = curve_to_q(X1);
%     [q2] = curve_to_q(X2);
%     [n] = size(q1,1);
%     A = q1*q2';
%     [U,S,V] = svd(A);
%     if det(A) > 0
%         S = eye(n);
%     else
%         S = eye(n);
%         S(:,end) = -S(:,end);
%     end
%     Ot = U*S*V';
%     X2 = Ot*X2;
%     q2 = Ot*q2;

% Applying optimal re-parameterization to the second curve
   
    