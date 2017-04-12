function [G] = snp_hhgeno(geno,ancalle)
%SNP_HHGENO - determins a SNP homozygote or heterozygote
% Syntax: [G] = snp_hhgeno(geno)
%
%disp('1   : Homozygote-Common allele ')
%disp('2   : Homozygote-Rare allele')
%disp('3   : Heterozygote')
%disp('4   : Undetermined ')

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-02-28 13:55:55 -0600 (Thu, 28 Feb 2013) $
% $LastChangedRevision: 462 $
% $LastChangedBy: jcai $

if nargin<2
    ancalle=[];
end

ddg=snp_ddgeno(geno);
ddg(ddg==55)=4;      	                 % Undetermined
ddg(~ismember(ddg,[11 22 33 44 4]))=3;   % Heterozygote

G=[];
mm=size(geno,2);
m=mm/2;

if ~isempty(ancalle)
    if length(ancalle)~=m
            error('Wrong ANCALLE provided.')
    end
end


if isempty(ancalle)
    for k=1:m
          x=ddg(:,k);
          y=x;
          y(y==3)=[];
          y(y==4)=[];
          [a,~,c]=unique(y);
          switch length(a)
              case 0
                          % do nothing if all individuals are heterozygotes
              case 1
                x(x==a(1))=1;       % Homozygote-Common (major) allele = 1
              case 2
                  if sum(c==1)>sum(c==2)
                      x(x==a(1))=1;
                      x(x==a(2))=2;   % Homozygote-Rare (minor) allele = 2
                  else
                      x(x==a(1))=2;
                      x(x==a(2))=1;
                  end                 
              otherwise
                  warning('SNP_HHGENO does not handle SNP with more than two alleles.')                  
                  x(:)=4;
          end
              G=cat(2,G,x);
    end
else
   % with ancalle
   disp('with anc')
    for k=1:m
      x=ddg(:,k);
      y=x;
      y(y==3)=[];
      y(y==4)=[];
      [a,~,c]=unique(y);
      ax=round(a./10);
      switch (length(a))
          case (0)
              %
          case (1)
              if ax(1)==ancalle(k)
                  x(x==a(1))=1;   % Homozygote-Common allele = 1
              else
                  x(x==a(1))=2;   % Homozygote-Common allele = 1
              end
          case (2)
              if ax(1)==ancalle(k)
                  x(x==a(1))=1;   % Homozygote-Common allele = 1
                  x(x==a(2))=2;   % Homozygote-Rare allele = 2
              elseif ax(2)==ancalle(k)
                  x(x==a(1))=2;
                  x(x==a(2))=1;
              else
                  x(x==a(1))=2;
                  x(x==a(2))=2;
                  warning('ANCALLE Cannot be determined.')
              end
          otherwise
%              error('SNP_HHGENO does not handle SNP with more than two alleles.')
      end
       	  G=cat(2,G,x);
    end
end
