function [nsite,ssite,isite]=geneid2nsynSITE(geneid)

nsite=-1;ssite=-1;isite=-1;
%transid='ENST00000360131';
%geneid='ENSG00000084693';
%genename='AGBL5';

[transidlist,transcord,D]=getcdscord(geneid);
%[transid2,transcord2,D]=genename2cds(genename);

if isempty(D)
    return;
end

    cdslen=zeros(1,length(transidlist));
    for (k=1:length(transidlist)),
          cdspos=transcord{k};
	      cdslen(k)=sum(cdspos(:,2)-cdspos(:,1));
    end
    [tempx,mxid]=max(cdslen);
    cdspos=transcord{mxid};
    
    
    chrid=D.chrid;
    strand=D.strand;
 	mmFilename=sprintf('C:/biodata/hgenome/Homo_sapiens.NCBI36.40.dna.chromosome.%d.mm',chrid);
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


[XS,XN]=getsynnonsynsites;

    ssite=sum(XS(codonise64(seq)));
    nsite=sum(XN(codonise64(seq)));
    
lentot=endn-startn+1;
lencds=sum(sum(cdspos(:,2)-cdspos(:,1)+1));
isite=lentot-lencds;



    