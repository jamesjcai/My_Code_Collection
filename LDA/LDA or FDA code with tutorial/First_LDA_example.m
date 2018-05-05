% prophet Mohammed said [ALLAH will help any one helped his/her brother/sister] PBUH

% This code is written by Alaa Tharwat (Frankfurt University of Applied Science - Germany)
% For more details about the code of the numerical example(s) see our paper "Tharwat, A., Gaber, T., Ibrahim, A., 
% & Hassanien, A. E. Linear discriminant analysis: A detailed tutorial. AI Communications,
% (Preprint), 1-22.?
% engalaatharwat@hotmail.com

clc
clear all
% This example deals with 2 classes
c1=[1 2;2 3;3 3;4 5;5 5]  % the first class (5 observations)
c2=[4 2;5 0;5 2;3 2;5 3;6 3] % the second class (6 observations)

% Number of observations of each class
n1=size(c1,1)
n2=size(c2,1)
N=n1+n2

%Mean of each class
mu1=mean(c1) % eq. (18)
mu2=mean(c2) % eq. (18)

% Average of the mean of all classes
mu=((n1/N)*mu1+(n2/N)*mu2) % eq. (18)

% Center the data (data-mean)
d1=c1-repmat(mu1,size(c1,1),1)
d2=c2-repmat(mu2,size(c2,1),1)


% Calculate the within class variance (SW)
s1=d1'*d1 % eq. (23)
s2=d2'*d2 % eq. (23)
sw=s1+s2 % eq. (23)
invsw=inv(sw) % eq. (24)

% in case of two classes only use v
% v=invsw*(mu1-mu2)'

% if more than 2 classes calculate between class variance (SB)
sb1=n1*(mu1-mu)'*(mu1-mu) % eq. (19)
sb2=n2*(mu2-mu)'*(mu2-mu) % eq. (20)
SB=sb1+sb2 % eq. (21)
v=invsw*SB % eq. (24)


% find eigne values and eigen vectors of the (v)
[evec,eval]=eig(v) % eq. (25)


% Sort eigen vectors according to eigen values (descending order) and
% neglect eigen vectors according to small eigen values
% v=evec(greater eigen value)
% or use all the eigen vectors

% project the data of the first and second class respectively
y1=c1*evec(:,2) % eq. (26)
y2=c2*evec(:,2) % eq. (27)

% You can check the other eigenvector (evec(:,1))
