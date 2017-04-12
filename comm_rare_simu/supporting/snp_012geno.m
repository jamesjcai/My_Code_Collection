function [G] = snp_012geno(geno)
%SNP_012GENO - Simplify GENO coding convention into 0, 1 and 2.

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[preG] = snp_hhgeno(geno);
%disp('1   : Homozygote-Common allele ')
%disp('2   : Homozygote-Rare allele')
%disp('3   : Heterozygote')
%disp('4   : Undetermined ')



%0 - AA   major allele
%1 - Aa
%2 - aa   minor allele
G=ones(size(preG));
G(preG==1)=0;   %disp('1   : Homozygote-Common allele ')
G(preG==2)=2;   %disp('2   : Homozygote-Rare allele')
G(preG==4)=nan;
