% prophet Mohammed said [ALLAH will help any one helped his/her brother/sister] PBUH
% This code to explain LDA or FDA with multiclass problems
% This code is written by Alaa Tharwat (Frankfurt University of Applied Science - Germany)
% For more details about the code of the numerical example(s) see our paper "Tharwat, A., Gaber, T., Ibrahim, A., 
% & Hassanien, A. E. Linear discriminant analysis: A detailed tutorial. AI Communications,
% (Preprint), 1-22.?
% engalaatharwat@hotmail.com
clc
clear all

% This example deals with 2 classes
c1=[0     1;
     0     2;
     1     2;
     2     0;
     1    1]  % the first class 5 observations

 c2=[ 1     9;
     2     9;
     1     7;
     2     8;
     1     8;] % the second class 6 observations

 c3=[5     8;
     7    9;
     6     10;
     5    7;
4	6] % the second class 6 observations

% Number of observations of each class
n1=size(c1,1)
n2=size(c2,1)
n3=size(c3,1)

%Mean of each class
mu1=mean(c1)
mu2=mean(c2)
mu3=mean(c3)


% Average of the mean of all classes
mu=(n1*mu1+n2*mu2+n3*mu3)/(n1+n2+n3)

% Center the data (data-mean)
d1=c1-repmat(mu1,size(c1,1),1)
d2=c2-repmat(mu2,size(c2,1),1)
d3=c3-repmat(mu3,size(c3,1),1)


% Calculate the within class variance (SW)
s1=d1'*d1
s2=d2'*d2
s3=d3'*d3
sw=s1+s2+s3
invsw=inv(sw)
% 
% % in case of two classes only use v
% v=invsw*(mu1-mu2)'

% if more than 2 classes calculate between class variance (SB)
sb1=n1*(mu1-mu)'*(mu1-mu)
sb2=n2*(mu2-mu)'*(mu2-mu)
sb3=n3*(mu3-mu)'*(mu3-mu)
SB=sb1+sb2+sb3
W=invsw*SB

% find eigne values and eigen vectors of the (v)
[evec,eval]=eig(W)

% Sort eigen vectors according to eigen values (descending order) and
% neglect eigen vectors according to small eigen values
% v=evec(greater eigen value)
% or use all the eigen vectors

% project the data of the first and second class respectively
y1=c1*evec(:,1)
y2=c2*evec(:,1)
y3=c3*evec(:,1)

% secoond eigenvector
y1_2=c1*evec(:,2)
y2_2=c2*evec(:,2)
y3_2=c3*evec(:,2)



%To plot PDF
% L1=min(y1):(max(y1)-min(y1))/20:max(y1)
% L2=min(y2):(max(y2)-min(y2))/20:max(y2)
% L3=min(y3):(max(y3)-min(y3))/20:max(y3)
% meany1=mean(y1)
% meany2=mean(y2)
% meany3=mean(y3)
% stdy1=std(y1)
% stdy2=std(y2)
% stdy3=std(y3)
% y1pdf=mvnpdf(L1',meany1,stdy1)
% y2pdf=mvnpdf(L2',meany2,stdy2)
% y3pdf=mvnpdf(L3',meany3,stdy3)
% histfit(min(y1):(max(y1)-min(y1))/20:max(y1),20),hold on
% histfit(min(y2):(max(y2)-min(y2))/20:max(y2),20), hold on
% histfit(min(y3):(max(y3)-min(y3))/20:max(y3),20)
% plot(y1(:,1),0,'r*','MarkerSize',10,'MarkerFaceColor','r')
% plot(y2(:,1),0,'bs','MarkerSize',10,'MarkerFaceColor','b')
% plot(y3(:,1),0,'go','MarkerSize',10,'MarkerFaceColor','g')
% legend('Class 1','Class 2','Class 3')



