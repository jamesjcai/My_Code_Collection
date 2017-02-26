function ci  = repli2ci(bstat,stat,alpha)
%ci  = repli2ci(bstat,stat,alpha)
%   CI = BOOTCI(NBOOT,BOOTFUN,...) computes the 95 percent BCa bootstrap
%   confidence interval of the statistic defined by the function BOOTFUN.
%   NBOOT is a positive integer indicating the number of bootstrap data
%   samples used in the computation. BOOTFUN is a function handle specified
%   with @. The third and later input arguments to BOOTCI are data
%   (scalars, column vectors, or matrices) that are used to create inputs
%   to BOOTFUN. BOOTCI creates each bootstrap sample by sampling with
%   replacement from the rows of the non-scalar data arguments (these must
%   have the same number of rows). Scalar data are passed to BOOTFUN
%   unchanged. CI is a vector containing the lower and upper bounds of the
%   confidence interval.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},'alpha',ALPHA) computes the
%   100*(1-ALPHA) percent BCa bootstrap confidence interval of the
%   statistic defined by the function BOOTFUN. ALPHA is a scalar between 0
%   and 1. The default value of ALPHA is 0.05.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'type',TYPE) computes the bootstrap
%   confidence interval of the statistic defined by the function BOOTFUN.
%   TYPE is the confidence interval type, specifying different methods of
%   computing the confidence interval. TYPE is a string chosen from
%       'norm' or 'normal':               normal approximated interval with
%                                         bootstrapped bias and standard
%                                         error;                                        
%       'per' or 'percentile':            basic percentile method; 
%       'cper' or 'corrected percentile': bias corrected percentile method;
%       'bca' :                           bias corrected and accelerated 
%                                         percentile method;
%       'stud' or 'student':              studentized confidence interval.
%   The default value of TYPE is 'bca'.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'type','stud','nbootstd',NBOOTSTD)
%   computes the studentized bootstrap confidence interval of the statistic
%   defined by the function BOOTFUN. The standard error of the bootstrap
%   statistics is estimated using bootstrap with NBOOTSTD bootstrap data
%   samples. NBOOTSTD is a positive integer value. The default value of
%   NBOOTSTD is 100.
%
%   CI = BOOTCI(NBOOT,{BOOTFUN,...},...,'type','stud','stderr',STDERR)
%   computes the studentized bootstrap confidence interval of statistics
%   defined by the function BOOTFUN. The standard error of the bootstrap
%   statistics is evaluated by the function STDERR. STDERR is a function
%   handle created using @. STDERR should take the same arguments as
%   BOOTFUN and return the standard error of the statistic computed by
%   BOOTFUN.
%
%   Example:
%     Compute the confidence interval for the capability index in
%     statistical process control:
%          y = normrnd(1,1,30,1);                  % simulated process data
%          LSL = -3;  USL = 3;                     % process specifications
%          capable = @(x) (USL-LSL)./(6* std(x));  % process capability
%          bootci(2000,capable, y)                 % Bca confidence interval
%          bootci(2000,{capable, y},'type','stud') % studentized confidence interval
%
%   See also: BOOTSTRP, JACKKNIFE.

% Copyright 2005 The MathWorks, Inc. 
 
    type = 'per';
if (nargin<3)
    alpha = .05;
end
 
% call sub-functions to compute the intervals 
switch (lower(type))
    case {'norm','normal'}
        ci = bootnorm(bstat,alpha,stat);
    case {'per','percentile'}
        ci = bootper(bstat,alpha);
    case {'cper', 'corrected percentile'}
        ci = bootcper(bstat,alpha,stat);        
    case 'bca'
        ci = bootbca(bstat,alpha,stat);
    otherwise 
        error('stats:bootci:BadType','BAD confidence interval type')
end;
 
%-------------------------------------------------------------------------    
function ci = bootnorm(bstat,alpha,stat)
% normal approximation interval
% A.C. Davison and D.V. Hinkley (1996), p198-200
 
se = std(bstat);   % standard deviation estimate
bias =mean(bstat-stat); % bias estimate
za = norminv(alpha/2);   % normal confidence point
lower = stat - bias + se*za; % lower bound
upper = stat - bias - se*za;  % upper bound
ci = [lower;upper];        
 
 
%-------------------------------------------------------------------------
function ci = bootper(bstat,alpha)
% percentile bootstrap CI 

pct1 = 100*alpha/2;
pct2 = 100-pct1;
lower = prctile(bstat,pct1); 
upper = prctile(bstat,pct2);
ci =[lower;upper];
 
%-------------------------------------------------------------------------
function ci = bootcper(bstat,alpha,stat)
% corrected percentile bootstrap CI
% B. Efron (1982), "The jackknife, the bootstrap and other resampling
% plans", SIAM.
 
% stat is transformed to a normal random variable z0.
% z0 = invnormCDF[ECDF(stat)]
x=length(bstat);
y=sum(bstat<stat);
z_0 = norminv(y/x);
z_alpha = norminv(alpha/2); % normal confidence point
 
% transform z0 back using the invECDF[normCDF(2z0-za)] and
% invECDF[normCDF(2z0+za)] 
pct1 = 100*normcdf(2*z_0-z_alpha); 
pct2 = 100*normcdf(2*z_0+z_alpha);
lower = prctile(bstat,pct2);  % inverse ECDF
upper = prctile(bstat,pct1);
ci = [lower;upper];
 
 
%-------------------------------------------------------------------------
function ci = bootbca(bstat,alpha,stat)
% corrected and accelerated percentile bootstrap CI
% T.J. DiCiccio and B. Efron (1996), "Bootstrap Confidence Intervals",
% statistical science, 11(3)
 
% same as bootcper, this is the bias correction
z_0 = norminv(sum(bstat<stat)./length(bstat));
 
% acceleration finding, see DiCiccio and Efron (1996)
jstat = jackknife(bootfun,varargin{:});
score = -(jstat-mean(jstat)); % score function at stat;
skew = sum(score.^3)./(sum(score.^2).^1.5);  %skewness of the score function
acc =  skew/6;  % accelleration
% tranform back with bias corrected and accelleration
z_alpha1 = norminv(alpha/2);
z_alpha2 = -z_alpha1;
pct1 = 100*normcdf(z_0 +(z_0+z_alpha1)/(1-acc*(z_0+z_alpha1)));
pct2 = 100*normcdf(z_0 +(z_0+z_alpha2)/(1-acc*(z_0+z_alpha2)));
% inverse of ECDF
ci = sort([prctile(bstat,pct2); prctile(bstat,pct1)]);
 