function [geno2] = snp_ddgeno(geno)
%SNP_DDGENO - double-digital encode genotype data.
% [1 3; 3 3] -> [13; 33]
% Syntax: [geno2] = snp_ddgeno(geno)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

mm=size(geno,2);
geno2=[];
for k=1:2:mm
    a=geno(:,k);
    b=geno(:,k+1);
    geno2=[geno2,a*10+b];
end