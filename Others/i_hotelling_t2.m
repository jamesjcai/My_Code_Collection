function [pvalue,T2]=i_hotelling_t2(X,mu)


[n,p]=size(X);
m=mean(X); %Mean vector from data matrix X.
S=cov(X);  %Covariance matrix from data matrix X.
%T2=n*(m-mu)*inv(S)*(m-mu)'; %Hotelling's T-Squared statistic.
T2=n*(m-mu)*(S\(m-mu)'); %Hotelling's T-Squared statistic.

if n >= 50 %Chi-square approximation.    
      X2=T2;
      v=p; %Degrees of freedom.
      pvalue=1-chi2cdf(X2,v); %Probability that null Ho: is true.
else  %F approximation.
      F=(n-p)/((n-1)*p)*T2;
      v1=p;  %Numerator degrees of freedom.
      v2=n-p;  %Denominator degrees of freedom.
      pvalue=1-fcdf(F,v1,v2);  %Probability that null Ho: is true.
end

% source: https://www.mathworks.com/matlabcentral/fileexchange/2844-hotellingt2?focused=5048346&tab=function
% https://onlinecourses.science.psu.edu/stat505/node/104
% https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc543.htm


