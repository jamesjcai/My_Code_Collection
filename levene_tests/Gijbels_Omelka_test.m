function [p,F,Dxv,Dyv]=Gijbels_Omelka_test(Dx,Dy)
% Testing for Homogeneity of Multivariate Dispersions Using
% Dissimilarity Measures
nx=size(Dx,1);
ny=size(Dy,1);

% Means of the winthin-group distances
dx_bar=mean(triu2vec(Dx));            % eq (2)
dy_bar=mean(triu2vec(Dy));            % eq (2)

% Average distance from observation to every other observation within the
% group
Dxv=sum(Dx,2)/(nx-1);             % eq (4)
Dyv=sum(Dy,2)/(ny-1);             % eq (4)

S2x=(4*(nx-1)/((nx-2)^2))*sum((Dxv-dx_bar).^2);     % eq (3)
S2y=(4*(ny-1)/((ny-2)^2))*sum((Dyv-dy_bar).^2);     % eq (3)
% sig_x=S2x/nx;  % eq (3)
% sig_y=S2x/ny;  % eq (3)

K=2;
n=nx+ny;
d_bar=(1/K)*(dx_bar+dy_bar);                 % eq (5)
sig=((nx-1)*S2x+(ny-1)*S2y)./(n-K);                % eq (5)

Fd=nx*(dx_bar-d_bar).^2+ny*(dy_bar-d_bar).^2;      % eq (5)
F=Fd./((K-1)*sig);                                 % eq (5)

df1=K-1; df2=n-K;
p=fpval(F,df1,df2);





