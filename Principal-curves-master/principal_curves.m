clear all
close all
% generating data
% dataflag == 1  spiral data simulation
% dataflag == 2  real data analysis
dataflag = 1;
if dataflag == 1
t = linspace(0,10*pi,500);
N1 = length(t); %sample size
n=N1;
x = 5*sin(t)+0.5*randn(n,1);
y = 5*cos(t)+0.5*randn(n,1);
z = t+0.5*randn(n,1);
data = [x,y,z]';
end
if dataflag == 2
    load sample.mat
    N1 = length(vol(1,1,:)); %sample size
    x = vol(:,1,1);
    y = vol(1,:,1);
    z = vol(1,1,:);
    data = [x,y,z]';
end

% PC projection with KDE pdf
data_in = data;
targetdim = 1; 
kernel_sigma = 1.5; 
pc_projection = pc_project_multidim(data,data_in,kernel_sigma,targetdim);
pc_projection = pc_projection';

figure, 
plot3(data(1,:),data(2,:),data(3,:),'.'); hold on,
title('PC projection with KDE pdf')
plot3(pc_projection(1,:),pc_projection(2,:),pc_projection(3,:),'.r'); axis equal;
