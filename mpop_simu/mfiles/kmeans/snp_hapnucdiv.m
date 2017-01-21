function [d]=snp_hapnucdiv(hap)
%SNP_HAPNUCDIV
%
%Syntax: [d]=snp_hapnucdiv(hap)

    [numHap,sizHap,seqHap]=counthaplotype(hap);
    p=sizHap./sum(sizHap);
    d=0;
    for i=1:numHap-1
    for j=i+1:numHap
        d=sum(seqHap(i,:)~=seqHap(j,:))./length(seqHap(j,:));
        % d=pdist(seqHap([i j],:),'hamming');
        d=d+2*p(i).*p(j).*d;
    end
    end
    d=d./(1-1/numHap);
