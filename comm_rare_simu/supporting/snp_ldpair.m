function [ldinfo] = snp_ldpair(genodata,markinfo,lenlimit)
%SNP_LDPAIR - calculate pairwise LD from genotype data
%
% [ldinfo] = snp_ldpair(genodata,markinfo,lenlimit)
% [ldinfo] = snp_ldpair(genodata)
%
% SNPs are from independent individuals
%
% EMLD is a program to calculate pair-wise linkage disequlibrium from
% genotype data. EM algorithm is used to estimate haplotype frequencies.
%
% SEE ALSO: EMLDRUN, SNP_LDPLOT

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% Email: jcai@tamu.edu
% 
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if (nargin<3), lenlimit=1000000; end
if (nargin<2), markinfo=[]; end

[n,m2]=size(genodata);
if (m2<4), error('Need more markers.'); end
m=m2/2;

ldinfo=struct;
ldinfo.d=zeros(m);
ldinfo.dprime=zeros(m);
ldinfo.r2=zeros(m);
warning off MATLAB:divideByZero;

if isempty(markinfo)
for i=1:m-1
for j=i+1:m
	    [d_raw,d_prime,r2,f]=i_calcld(i_getsnpgeno([i j],genodata));
	    ldinfo.d(i,j)=d_raw;
	    ldinfo.dprime(i,j)=abs(d_prime);
	    ldinfo.r2(i,j)=r2;
        ldinfo.freq{i,j}=f;
end
end
else
    for i=1:m-1
    for j=i+1:m
        if ((abs(markinfo.pos(i)-markinfo.pos(j)))<lenlimit),
            [d_raw,d_prime,r2,f]=i_calcld(i_getsnpgeno([i j],genodata));
            ldinfo.d(i,j)=d_raw;
            ldinfo.dprime(i,j)=abs(d_prime);
            ldinfo.r2(i,j)=r2;
            ldinfo.freq{i,j}=f;
        end
    end
    end
end
warning on MATLAB:divideByZero




%%%%%%%%%%%%
%%% SUB  %%%
%%%%%%%%%%%%

function [d_raw,d_prime,r2,probHaps] = i_calcld(snppair)
%Internal function

%    * L1 and L2 are the two loci in question, referenced by their number or name (if marker info file is provided)
%    * D' is the value of D prime between the two loci.
%    * LOD is the log of the likelihood odds ratio, a measure of confidence in the value of D'
%    * r2 is the correlation coefficient between the two loci
%    * CIlow is 95% confidence lower bound on D'
%    * CIhi is the 95% confidence upper bound on D'
%    * Dist is the distance (in bases) between the loci, and is only displayed if a marker info file has been loaded
%    * T-int is a statistic used by the HapMap Project to measure the completeness of information represented by a set of markers in a region

[probHaps,probMajor]=snp_hapfreqem(snppair);

fA=probMajor(1); fa=1-fA;
fB=probMajor(2); fb=1-fB;


%d_raw1=sum(diag(probHaps))./sum(probHaps(:))-prod(probMajor);
d_raw=probHaps(1,1)-prod(probMajor);

r2=(d_raw*d_raw)./(prod(probMajor)*prod(1-probMajor));
r2=min([1 r2]);

		if d_raw<0
            x=min([prod(probMajor),prod(1-probMajor)]);
        else
            x=min([fA*fb,fa*fB]);
        end        
        if x==0
            d_prime=nan;
            r2=nan;
        else
	        d_prime=min([1 d_raw./x]);
        end



function [sgeno]=i_getsnpgeno(index, geno, genoinfo)
%GETSNPGENO - extracts the index-th SNP's genotype data
%for example, i_getsnpgeno(5,geno)==geno(:,[9,10])

m=size(geno,2);
if (index>m/2), error ('Geno does not have so many SNP.'); end

sgeno=geno(:,i_newindex(index));

if (nargin>2),
    genoinfo.rsid([1 2])
end

function [ii]=i_newindex(i)
    %converts [1 2] => [1 2 3 4]
    m=length(i);
    ii=[];
    x=i*2-1;
    y=i*2;
    for k=1:m
      ii=[ii,x(k),y(k)];
    end
