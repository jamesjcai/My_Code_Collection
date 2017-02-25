function [nn,ns,nsynsite,synsite]=subs_at_ensemblgene(geneid)

nsynsite=-1;
synsite=-1;
nn=-1; ns=-1;
%transid='ENST00000360131'; %geneid='ENSG00000084693'; %genename='AGBL5';

[transidlist,transcord,D]=getcdscord(geneid);
%[transid2,transcord2,D]=genename2cds(genename);


if isempty(D)
    return;
end
[cdspos]=i_getlongestcds(transidlist,transcord);
[seqHs]=i_getseq(cdspos,D,'human');
[seqPt]=i_getseq(cdspos,D,'chimp');

seqHs(seqHs>4)=5;
seqPt(seqPt>4)=5;

seqori=[seqHs;seqPt];
seq=rmcodongaps(seqori);


if size(seq,2)/size(seqori,2)<0.7
    warning('More than 30% of orginal CDS was removed.')
end


if ~isempty(seq)
    [dS,dN,dN_dS,VdS,VdN,synsite,nsynsite,ns,nn]=dc_ng86(seq);
    ns=ns(1,2);
    nn=nn(1,2);
    synsite=synsite(1,2);
    nsynsite=nsynsite(1,2);
end





function [cdspos]=i_getlongestcds(transidlist,transcord)
    cdslen=zeros(1,length(transidlist));
    for (k=1:length(transidlist)),
          cdspos=transcord{k};
	      cdslen(k)=sum(cdspos(:,2)-cdspos(:,1));
    end
    [tempx,mxid]=max(cdslen);
    cdspos=transcord{mxid};
    
    
function [seq]=i_getseq(cdspos,D,species)
    chrid=D.chrid;
    strand=D.strand;
    
    switch species
        case {'human'}
         	mmFilename=sprintf('C:/biodata/hgenome/Homo_sapiens.NCBI36.40.dna.chromosome.%d.mm',chrid);
        case {'chimp'}
            mmFilename=sprintf('Y:/chimpaln/mmFilenamechr%d',chrid);
    end
	chrmm = memmapfile(mmFilename, 'format', 'uint8');    
    seq=[];
    for (kk=1:size(cdspos,1)),
          seq=[seq,chrmm.Data(cdspos(kk,1):cdspos(kk,2))'];
    end
    if (strand<0)
        seq=revcomseq(seq);
    end

startn=min(min(cdspos));
endn=max(max(cdspos));
