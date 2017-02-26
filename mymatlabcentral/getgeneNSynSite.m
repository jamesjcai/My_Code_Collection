function [ssitex,nsitex,transid,transcord,cdsseq]=getgeneNSynSite(geneid,chrid,strand,chrmm)

if (nargin<4)
mmFilename=sprintf('C:/biodata/hgenome/Homo_sapiens.NCBI36.40.dna.chromosome.%d.mm',chrid);
chrmm = memmapfile(mmFilename, 'format', 'uint8');
end


[transid,transcord]=getcdscord(geneid);
[XS,XN]=getsynnonsynsites;



nsitex={};
ssitex={};
cdsseq={};
for ss=1:length(transid)
    cdspos=transcord{ss};
    seq=[];
    for (k=1:size(cdspos,1)),
       seq=[seq,chrmm.Data(cdspos(k,1):cdspos(k,2))'];
    end
    if (strand<0)
	seq=revcomseq(seq);
    end

    zseq=codonise64(seq);
    if (any(zseq>64))
        idxz=find(zseq>64);
        for st=1:length(idxz)
            zseq(idxz(st))=floor(rand*64)+1;
        end
    end
    
    ssite=XS(zseq);
    nsite=XN(zseq);

    ssitex{ss}=ssite;
    nsitex{ss}=nsite;
    cdsseq{ss}=zseq;    
end
%clear chrmm;


    