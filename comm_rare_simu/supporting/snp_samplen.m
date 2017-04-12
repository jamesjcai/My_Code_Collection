function [smpln,markn,indvn,smplnv]=snp_samplen(geno)
%SNP_SAMPLEN - number of samples (chromosomes) in given genotype data
%
% [smpln,markn,indvn,smplnv]=snp_samplen(geno)
%
%See also: SNP_MARKLEN, SNPMARKBPLEN

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

[n,m2]=size(geno);
markn=m2/2;
smpln=n*2;     % number of chromosomes
indvn=n;       % number of individual



if nargin>3
    smplnv=zeros(1,markn);
    G=[geno(1:2:m2);geno(2:2:m2)];
    for k=1:markn
         smplnv(k)=sum(G(:,k)~=5);
    end
end