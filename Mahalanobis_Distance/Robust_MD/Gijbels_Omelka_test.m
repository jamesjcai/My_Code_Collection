function [p,F]=Gijbels_Omelka_test(Dx,Dy)
% Testing for Homogeneity of Multivariate Dispersions Using
% Dissimilarity Measures
nx=size(Dx,1);
ny=size(Dy,1);
dx_bar=mean(triu2num(Dx));
dy_bar=mean(triu2num(Dy));

Dxv=mean(Dx)*(nx/(nx-1));
S2x=(4*(nx-1)/(nx-2)^2)*sum((Dxv-dx_bar).^2);

Dyv=mean(Dy)*(ny/(ny-1));
S2y=(4*(ny-1)/(ny-2)^2)*sum((Dyv-dy_bar).^2);

k=2;

d_bar=0.5*(nx*dx_bar+ny*dy_bar);
sig=((nx-1)*S2x+(ny-1)*S2y)./(nx+ny-k);

Fd=nx*(dx_bar-d_bar).^2+ny*(dy_bar-d_bar).^2;
F=Fd./(k-1)*sig;

df1=k-1; df2=nx+ny-k;
p = fpval(F,df1,df2);





