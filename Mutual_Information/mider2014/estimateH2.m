%--------------------------------------------------------------------------
% estimateH2: obtain joint entropy and mutual information of two variables
% 
%   [mutinfo,fracn,H2] = estimateH2(x,y,pb,q) calculates the joint entropy
%   'H2' and mutual information 'mutinfo' of two variables 'x' and 'y', 
%   using the type of entropy specified by parameter 'q', and estimating the 
%   entropy with 'pb' data points per bin. It also checks if criteria for
%   X^2 (chi-square) test are met, and calculates the X^2 probability to 
%   reject null hypothesis.
% 
%   x       = data vector of the 1st variable
%   y       = data vector of the 2nd variable 
%   pb      = estimated number of points per bin
%   q       = entropic parameter; if q>1 it is Tsallis entropy
%   mutinfo = average mutual information of (x,y)
%   fracn   = fraction of joint prob cells that contain >= 5 points.
%   H2      = joint entropy of (x,y)
% 
%   Dependencies: estimateH2.m is called by mider.m
%
%	Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
%	Created: 07/10/2011, last modified: 17/04/2013
%--------------------------------------------------------------------------
% Copyright (C) 2013  Alejandro Fernández Villaverde
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------


function [mutinfo,fracn,H2] = estimateH2(x,y,pb,q)

szx     = size(x); if szx(1)>szx(2), x=x'; end  
szy     = size(y); if szy(1)>szy(2), y=y'; end
leng    = length(x);
binsize = ceil(sqrt(pb*leng));   % Estimated bin size for ~pb points/bin
b       = floor(leng/binsize);   % Number of bins
binsize = floor(leng/b);         % Integer number of elements per bin

% Initialize matrices to zero
his = zeros(b,b);
x1  = zeros(b,1);
x2  = zeros(b,1);

% Sort values of X and Y (descending order)
[sx,ix] = sort(x,2,'descend');   
[sy,iy] = sort(y,2,'descend');  

% Replace elements of X and Y by their rank order
linsp  = linspace(1,leng,leng);
sx(ix) = linsp;
sy(iy) = linsp;

% Put sorted elements into bins
sx = ceil(sx/binsize);
sy = ceil(sy/binsize);

% Put largest values in highest bin
for p=1:leng
   if sx(p)>b
        sx(p)=b;
   end
   if sy(p)>b
        sy(p)=b;
   end
end

% Populate histograms
for i=1:leng,
    x = sx(i);
    y = sy(i);
    his(y,x) = (his(y,x)+1);
    x1(x) = (x1(x)+1);
    x2(y) = (x2(y)+1);
end

% Calculate Expected Values
expct  = (x2*x1')/(sum(sum(his)));
dif    = expct-his;
chisqr = (dif.*dif)./expct;
X2     = sum(sum(chisqr)); % The chi-square value
df     = (b-1)^2; % Degrees of freedom used for the X^2 calculations
                     
% Continue with calculation of Mutual Information and H2
his = his/leng;
x1  = x1/leng;
x2  = x2/leng;
res = 0; 
res2 = 0;
if q>1
    for j = 1:b,
        for k = 1:b
            if(his(j,k)>0)
                an2  = his(j,k)*(his(j,k)^(1-q)-1)/(1-q);               % H2
                res2 = res2 + an2;                                      % H2
                an   = his(j,k)*(his(j,k)/(x1(k).*x2(j))^(1-q)-1)/(1-q);% MI
                res  = res + an;                                        % MI 
            end
        end
    end
else
    for j = 1:b,
        for k = 1:b
            if(his(j,k)>0)
                an2  = his(j,k)*log2(his(j,k));                         % H2
                res2 = res2 + an2;                                      % H2
                an   = his(j,k)*log2(his(j,k)/(x1(k).*x2(j)));          % MI
                res  = res + an;                                        % MI
            end 
        end
    end
end
                                        
mutinfo = res;
H2      = -res2;

% Check to see if the chi-square criteria are met (chi-square probability)
prob = chi2pdf(X2,df); 
    % prob = 1 minus the probability that a single observation from a 
    % chi square distribution with df degrees of freedom will fall within 
    % the interval [0 X2] (e.g. if prob < 0.01, reject the null hypothesis)
% If the hypothesis of statistical independence is not rejected, make I(X,Y)=0
if prob > 0.01, mutinfo = 0; end 

% Calculate fracn
[a,b]   = size(his);
his     = his*leng;
m =1 ; 
for k=1:a;
    for j=1:b;
        if his(k,j)>0,
            hnz(m)=his(k,j);
            m=m+1;
        end;
    end;
end
lenz = length(hnz); 
num5 = 0; for k=1:lenz;if hnz(k)>=5;num5=num5+1;end;end
fracn = num5/lenz;
