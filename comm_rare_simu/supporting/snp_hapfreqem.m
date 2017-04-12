function [probHaps,probMajor,seqHap1,seqHap2]=snp_hapfreqem(snppair)
%SNP_HAPFREQEM - EM algorithm to estimate probabilities of haplotypes
% Syntax: [probHaps,probMajor]=snp_hapfreqem(snppair)

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $


[~,m]=size(snppair);
if (m~=4), error('Need SNP pair, m=4'); end
snppair(sum(snppair==5,2)>0,:)=[];   % remove genotype=5


%[numHap,sizHap,seqHap] = counthaplotype(snppair);
[x3,seqHap1,seqHap2]=i_build3x3table(snppair);
[~,~,probHaps]=i_3to2(x3);

%[probHaps] = probhaps(x2,xhet)
%[probMajor]=1-snp_maf(snppair)
[p1,~,q1]=i_pq(x3,size(snppair,1));
probMajor=[p1,q1];

% loglike = i_loglike below
%        loglike  = (known[AA]*Math.log(probHaps[AA]) + known[AB]*Math.log(probHaps[AB]) + known[BA]*Math.log(probHaps[BA]) + known[BB]*Math.log(probHaps[BB]))/LN10 + ((double)unknownDH*Math.log(probHaps[AA]*probHaps[BB] + probHaps[AB]*probHaps[BA]))/LN10;
%        loglike1 = (known[AA]*Math.log(probHaps[AA]) + known[AB]*Math.log(probHaps[AB]) + known[BA]*Math.log(probHaps[BA]) + known[BB]*Math.log(probHaps[BB]) + (double)unknownDH*Math.log(probHaps[AA]*probHaps[BB] + probHaps[AB]*probHaps[BA]))/LN10;
%        loglike0 = (known[AA]*Math.log(pA1*pA2) + known[AB]*Math.log(pA1*pB2) + known[BA]*Math.log(pB1*pA2) + known[BB]*Math.log(pB1*pB2) + (double)unknownDH*Math.log(2*pA1*pA2*pB1*pB2))/LN10;
% LOD=loglike1-loglike0;
% ref: HaploData.java
% D' is the value of D prime between the two loci.
% LOD is the log of the likelihood odds ratio, a measure of confidence in
% the value of D'.


function [x3,seqHap1,seqHap2]=i_build3x3table(snppair)
    x3=zeros(3);
    n=size(snppair,1);
    %      AA  |  Aa  |  aa
    %    -------------------
    % BB |  0  |   1  |  2
    % Bb |  3  |   4  |  5
    % bb |  6  |   7  |  8

    % A = major allele; B = major allele
    [s1,seqHap1]=i_encodesnp(snppair(:,[1 2]),n);
    [s2,seqHap2]=i_encodesnp(snppair(:,[3 4]),n);
    % AA=BB=1; Aa=Bb=2; aa=bb=3
    for k=1:n
        x3(s1(k),s2(k))=x3(s1(k),s2(k))+1;
    end



function [snpx,seqHap]=i_encodesnp(s1,n)
    [numHap,~,seqHap] = counthaplotype(s1(:));
    if (numHap>2), error('Something wrong'); end
    snpx=2*ones(n,1);
    for k=1:n
        if (s1(k,1)==seqHap(1) && s1(k,2)==seqHap(1)),
            snpx(k)=1;
        elseif (s1(k,1)==seqHap(2) && s1(k,2)==seqHap(2)),
            snpx(k)=3;
        end
    end


function [p1,p2,q1,q2]=i_pq(x3,n)
		p1 = (2*(x3(1,1)+x3(1,2)+x3(1,3))+x3(2,1)+x3(2,2)+x3(2,3))/(2*n);
		p2 = (2*(x3(3,1)+x3(3,2)+x3(3,3))+x3(2,1)+x3(2,2)+x3(2,3))/(2*n);
		q1 = (2*(x3(1,1)+x3(2,1)+x3(3,1))+x3(1,2)+x3(2,2)+x3(3,2))/(2*n);
        if nargout>3
		q2 = (2*(x3(1,3)+x3(2,3)+x3(3,3))+x3(1,2)+x3(2,2)+x3(3,2))/(2*n);
        end



function [x2,xhet,probHaps]=i_3to2(x3)

%===============
%  |  A  |  a  |
%   ------------
%B | 1,1 | 1,2 |
%b | 2,1 | 2,2 |
%===============
%  nAB=(float)(2*cells[0]+cells[1]+cells[3]);
%  nab=(float)(2*cells[8]+cells[7]+cells[5]);
%  nAb=(float)(2*cells[6]+cells[7]+cells[3]);
%  naB=(float)(2*cells[2]+cells[1]+cells[5]);
%  nAB = 2*AABB + AaBB + AABb;
%  nab = 2*aabb + Aabb + aaBb;
%  nAb = 2*AAbb + Aabb + AABb;
%  naB = 2*aaBB + AaBB + aaBb;

x2=zeros(2,2);
x2(1,1)=2*x3(1,1)+x3(1,2)+x3(2,1);   % AB
x2(2,2)=2*x3(3,3)+x3(2,3)+x3(3,2);   % ab
x2(1,2)=2*x3(1,3)+x3(2,3)+x3(1,2);   % Ab
x2(2,1)=2*x3(3,1)+x3(2,1)+x3(3,2);   % Ba

% 1/22/2010: a bug at above two lines has been fixed.
% thanks Eran Elhaik for pointing it out.



xhet=x3(2,2);
N=sum(x2(:))+2*xhet;
probHaps=zeros(2);
dx=100;
theta=-999999999.0;
count=1;
while (count<1000 && dx>1e-8),
    thetaprev=theta;
    count=count+1;
    prAB=(x2(1,1) + (1-theta)*xhet)/N;
    prab=(x2(2,2) + (1-theta)*xhet)/N;
    prAb=(x2(2,1) + theta*xhet)/N;
    praB=(x2(1,2) + theta*xhet)/N;
    theta=(prAb*praB)/(prAB*prab + prAb*praB);
    dx=abs(theta-thetaprev);
end
probHaps(1,1)=prAB;
probHaps(2,2)=prab;
probHaps(2,1)=prAb;
probHaps(1,2)=praB;
%assert(2*sum(x3(:))==(2*xhet + sum(x2(:))))



%{
function [probHaps] = probhaps(knownHaps,unknownDH)
%EM algorithm to estimate probabilities of haplotypes
%EM Computation of Haplotype Probabilities

%knownHaps=[122,0;0 18];
%knownHaps=[191 0;3 30];

numHaps=knownHaps;
const_prob=0.1;
probHaps=i_estimateP(knownHaps,const_prob);

% now we have an initial reasonable guess at p we can
% start the EM

%const_prob=0.0;
count=1; loglike=-999999999.0;
TOLERANCE=0.00000001;
%TOLERANCE=0.01;

%unknownDH=4;

x=10;
while ((count<1000)&&(x > TOLERANCE)),
	oldloglike=loglike;
	loglike=i_loglike(knownHaps,probHaps,unknownDH);
	x=abs(loglike-oldloglike);
	numHaps2=i_count_haps(numHaps,probHaps,count,knownHaps,unknownDH);
	numHaps=numHaps2;
	probHaps=i_estimateP(numHaps,0);
	count=count+1;
end
%disp(count);
%}



%~*~*~*~*~*~*~*~*~*
%  SUBFUNCTIONS   *
%~*~*~*~*~*~*~*~*~*

%{
function [like] = i_loglike(knownHaps,probHaps,unknownDH)
	like=sum(sum(knownHaps.*log(probHaps))) +...
		(unknownDH*log(prod(diag(probHaps)) + probHaps(1,2)*probHaps(2,1)));


function [P] = i_estimateP(N,const_prob)
	total=sum(N(:))+4.0*const_prob;
	P=(N+const_prob)./total;
	P(P<1e-10)=1e-10;


function [numHaps2] = i_count_haps(numHaps,probHaps,em_round,knownHaps,unknownDH)
	numHaps2=knownHaps;
	if (em_round > 0),
            numHaps2(1,1) = numHaps2(1,1) + unknownDH* (probHaps(1,1)*probHaps(2,2))/((probHaps(1,1)*probHaps(2,2))+(probHaps(1,2)*probHaps(2,1)));
            numHaps2(2,2) = numHaps2(2,2) + unknownDH* (probHaps(1,1)*probHaps(2,2))/((probHaps(1,1)*probHaps(2,2))+(probHaps(1,2)*probHaps(2,1)));
            numHaps2(1,2) = numHaps2(1,2) + unknownDH* (probHaps(1,2)*probHaps(2,1))/((probHaps(1,1)*probHaps(2,2))+(probHaps(1,2)*probHaps(2,1)));
            numHaps2(2,1) = numHaps2(2,1) + unknownDH* (probHaps(1,2)*probHaps(2,1))/((probHaps(1,1)*probHaps(2,2))+(probHaps(1,2)*probHaps(2,1)));
	end
%}