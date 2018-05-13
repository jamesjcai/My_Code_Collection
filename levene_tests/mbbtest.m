function mbbtest(X,t,alpha)
%MBBTEST Multivariate Bootstrap Bartlett´s Test for the Homogeneity of 
%Covariance Matrices.
%  The bootstrap is a way of estimating the variability of a statistic   
%  from a single data set by resampling it independently and with equal
%  probabilities (Monte Carlo resampling). Allows the estimation of 
%  measures where the underlying distribution is unknown or where sample 
%  sizes are small. Their results are consistent with the statistical 
%  properties of those analytical methods (Efron and Tibshirani, 1993).
%
%  The name 'bootstrap' originates from the expression 'pulling yourself 
%  up by your own bootstraps' and refers to the basic idea of the 
%  bootstrap, sampling with replacement from the data. In this way a
%  large number of 'bootstrap samples' is generated, each of the same size
%  as the original data set. From each bootstrap sample the statistical 
%  parameter of interest is calculated (Wehrens and Van der Linden, 1997).
%  
%  Here, we use the Non-parametric Bootstrap. Non-parametric bootstrap is 
%  simpler. It does not use the structure of the model to construct 
%  artificial data. The data is instead directly resampled with
%  replecement.
%
%  As Silva et. al (2008) stated, here a m-file analytical procedure using
%  bootstrap method is developed as an alternative to the homogeinity of
%  covariance matrices test. 
%
%  MBBTEST treats NaN values as missing values, and removes them.
%
%  Syntax: function mbbtest(X,t,alpha) 
%      
%  Inputs:
%       X - data matrix (Size of matrix must be n-by-(1+p); sample=
%           column 1, variables=column 2:p)
%       t - boot times or number of Bootstrap simulations (resamplings)
%   alpha - significance level (default = 0.05)
%  Output:
%         - Whether or not the homoscedasticity was met
%  
%  Example: We want to assess the homoscedasticity of covariance matrices
%           of the effects of cognitive behaviour therapy (CBT=G1) on 
%           obsessive compulsive disorder (OCD). CBT will we compared with
%           Behavior Therapy (BT=G2) and with no-treatment (NT=G3) as a  
%           control condition. Since OCD manifests itself both behaviorally 
%           (obsessive actions=x1) as well as cognitively (obsessive 
%           thoughts=x2), both will be measured. Note that the two dependent
%           variables are theoretically motivated. We have three groups and
%           two dependent variables. It was taken from online at:
%           http://www.ii.metu.edu.tr/~hohenberger/statistics/Chapter_14.ppt
% 
%                                    Group
%                    ---------------------------------------                
%                          1           2           3
%                    ---------------------------------------
%                       x1    x2    x1    x2    x1    x2
%                    ---------------------------------------
%                        5    14     4    14     4    13
%                        5    11     4    15     5    15
%                        4    16     1    13     5    14
%                        4    13     1    14     4    14
%                        5    12     4    15     6    13
%                        3    14     6    19     4    20
%                        7    12     5    13     7    13
%                        6    15     5    18     4    16
%                        6    16     2    14     6    14
%                        4    11     5    17     5    18
%                    ---------------------------------------
%
%  Total data matrix must be:
%  X=[1 5 14;1 5 11;1 4 16;1 4 13;1 5 12;1 3 14;1 7 12;1 6 15;1 6 16;1 4 11;
%     2 4 14;2 4 15;2 1 13;2 1 14;2 4 15;2 6 19;2 5 13;2 5 18;2 2 14;2 5 17;
%     3 4 13;3 5 15;3 5 14;3 4 14;3 6 13;3 4 20;3 7 13;3 4 16;3 6 14;3 5 18];
%
%  Calling on Matlab the function: 
%            mbbtest(X,1000,0.05)
%
%  Answer is:
%
%  The number of groups are: 3
%  The number of dependent variables are: 2
%
%  Bootstraping probability associated = 0.1810
%  The bootstrap associated probability is equal or larger than 0.05
%  So, after 1000 bootstraps the assumption of homoscedasticity was met.
%
%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.edu.mx
%
%  Copyright (C)  November 30, 2011
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2011). mbbtest:Multivariate
%     Bootstrap Bartlett´s Test for the Homogeneity of Covariance Matrices.
%     [WWW document].URL http://www.mathworks.com/matlabcentral/fileexchange/
%     34037-mbbtest
%
%  References:
%  Efron, B. and Tibshirani, R. J. (1993), An Introduction to the Bootstrap
%             Chapman and Hall:New York.
%  Silva, R. B. V., Ferreira, D. F. and Nogueira, D. A. (2008), Robustness 
%             of asymptotic and bootstrap tests for multivariate homogeneity
%             of covariance matrices. Ciênc. agrotec., Lavras, 32(1):157-166.
%             Available online at:http://www.scielo.br/scielo.php?pid=S1413-
%             70542008000100023&script=sci_arttext
%  Wehrens, R and Van der Linden, W. E. (1997), Bootstrapping Principal
%             Component Regression Models. Journal of Chemometrics, 
%             11:157–171.
%

if  nargin < 2,
    error('mbbtest:TooFewInputs', ...
          'MBBTEST requires at least two input arguments.');
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05; %default
elseif numel(alpha) ~= 1 || alpha <= 0 || alpha >= 1
    error('mbbtest:BadAlpha','ALPHA must be a scalar between 0 and 1.');
end

%Remove NaN values, if any
X = X(~any(isnan(X),2),:);

g = max(X(:,1)); %Number of groups
fprintf('The number of groups are:%2i\n', g);
c = size(X,2);
p = c-1; %Number of dependent variables
fprintf('The number of dependent variables are:%2i\n', p);
disp(' ')

xx = X(:,1);

indice = X(:,1);
for i = 1:g
    Xe = indice == i;
    s(i).X = X(Xe,2:c);
    s(i).m = mean(s(i).X);
    s(i).cv = cov(s(i).X);
    s(i).n = length(s(i).X);
    s(i).df = s(i).n - 1;
    s(i).a = s(i).df*s(i).cv;
    s(i).b = s(i).df*log(abs(det(s(i).cv))); %abs to avoid a complex number
    s(i).e = s(i).X - repmat(s(i).m,s(i).n,1);
end
n=cat(1,s.n);df=cat(1,s.df);b=cat(1,s.b);e=cat(1,s.e);
%e=pooled residuals (the order does not matter)

N = sum(n);

sm = zeros(size(s(1).a));
for i = 1:g
    sm = sm + s(i).a;
end

Sp = sm/(N-g);
C = (N - g)*log(det(Sp)) - sum(b); %Box's M statistic
X2 = (1-(sum(1./df)-(1/(N-g)))*(((2*p^2)+(3*p)-1)/(6*(p+1)+(g-1))))*C;

idx = ceil(rand(N,t)*N);

TB = [];
for i = 1:t
    X = [xx e(idx(:,i),:)];
    for j = 1:g
        Xe = indice == j;
        s(j).X = X(Xe,2:c);
        s(j).m = mean(s(j).X);
        s(j).cv = cov(s(j).X);
        s(j).n = length(s(j).X);
        s(j).df = s(j).n - 1;
        s(j).a = s(j).df*s(j).cv;
        s(j).b = s(j).df*log(abs(det(s(j).cv))); %abs to avoid a complex 
    end                                           %number
    n=cat(1,s.n);df=cat(1,s.df);b=cat(1,s.b);
    
    N = sum(n);
    
    sm = zeros(size(s(1).a));
    for j = 1:g
        sm = sm + s(j).a;
    end
    
    Sp = sm/(N-g);
    C = (N - g)*log(det(Sp)) - sum(b); %Box's M statistic
    X2b = (1-(sum(1./df)-(1/(N-g)))*(((2*p^2)+(3*p)-1)/(6*(p+1)+(g-1))))*C;
    TB = [TB;X2b];
end

X2c = length(find(TB > X2));
PX2 = X2c/t; %P-value associated to the bootstrap Chi-squared statistic

fprintf('Bootstraping probability associated = %3.4f\n', PX2);

if PX2 >= alpha;
  fprintf('The bootstrap associated probability is equal or larger than% 3.2f\n', alpha);
  fprintf('So, after %i bootstraps the assumption of homoscedasticity was met.\n', t);
else
  fprintf('The bootstrap associated probability is smaller than% 3.2f\n', alpha);
  fprintf('So, after %i bootstraps the assumption of homoscedasticity was not met.\n', t);
end

return,