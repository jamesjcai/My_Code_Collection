function [numHap,sizHap,seqHap,idxHap] = counthaplotype(aln)
%COUNTHAPLOTYPE - Counts haplotypes

% Population Genetics and Evolution Toolbox (PGEToolbox)
% Author: James Cai
% (c) Texas A&M University
%
% $LastChangedDate: 2013-01-06 13:39:38 -0600 (Sun, 06 Jan 2013) $
% $LastChangedRevision: 331 $
% $LastChangedBy: jcai $

if isstruct(aln), seq=aln.seq; else seq=aln; end
[n,m]=size(seq);
if n==1, numHap=1; sizHap=1; seqHap=seq; return; end

if nargout>3
    methodid=1;
else
    methodid=2;
end

switch methodid
    case 1
        [x,y,z]=unique(seq,'rows');
        numHap=size(x,1);

        if nargout>1
            sizHap=grpstats(z,z,'numel');
            [sizHap,y2]=sort(sizHap,'descend');
            %sizHap=sizHap(end:-1:1);
        end
        if nargout>2
            seqHap=seq(y,:);
            %seqHap=seqHap(end:-1:1,:);
            seqHap=seqHap(y2,:);
        end
        if nargout>3
            
           %idxHap=-1*(z-(max(z)+1));
           idxHap=zeros(size(z));
           for k=1:numHap
               idxHap(z==y2(k))=k;
           end
        end
        

    case 2
        % faster
        [seq,idx]=sortrows(seq);
        seqHap=seq(1,:);
        curseq=seq(1,:);
        sizHap=zeros(n,1);
        numHap=1;
        for i=1:n
            if  sum(seq(i,:)==curseq)~=m
            	seqHap=[seqHap;seq(i,:)];
                numHap=numHap+1;
            	curseq=seq(i,:);
            end
            sizHap(numHap)=sizHap(numHap)+1;
        end
        if nargout>1
          sizHap(numHap+1:end)=[];
          [sizHap,idx]=sort(sizHap,'descend');
        end
        if nargout>2
           seqHap=seqHap(idx,:);
        end
end


