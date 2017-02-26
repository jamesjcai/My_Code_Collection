function [nn,ns,ni,D]=dbsnp128_at_ensemblexon(geneid,transid,exonstart,exonend)

nn=-1; ns=-1; ni=-1;
%transid='ENST00000360131';
%geneid='ENSG00000084693';
%genename='AGBL5';
%[nn,ns,ni,D]=dbsnp128_at_ensemblexon('ENSG00000188157','ENST00000379364',976275,976286)
%[nn,ns,ni,D]=dbsnp128_at_ensemblexon('ENSG00000116731','ENST00000235372',13977500,13981913)

[transidlist,transcord,D]=getcdscord(geneid);
%[transid2,transcord2,D]=genename2cds(genename);

if isempty(D)
    return;
end

[ok,mxid]=ismember(transid,transidlist);
if ~ok
    return;
end
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


%startn=min(min(cdspos));
%endn=max(max(cdspos));

%startn=min(min(cdspos));
%endn=max(max(cdspos));
startn=exonstart;
endn=exonend;

    % get SNPs
  
%sprintf(' WHERE (((r.chr)=''%s'') AND ((r.chrpos)>=%d And (r.chrpos)<=%d) AND ((r.totalhits)<2));'...

sql='SELECT r.* FROM Refchr_all r';
sql=strcat(sql,...
sprintf(' WHERE (((r.chr)=''%s'') AND ((r.chrpos)>=%d AND (r.chrpos)<=%d));'...
,num2str(chrid),startn,endn));


conn = database('dbsnp128', '', '');
curs = exec(conn, sql);
setdbprefs('DataReturnFormat','structure')
curs = fetch(curs);
D=curs.Data;
close(curs);
close(conn);
nn=0; ns=0; ni=0;
if ~isstruct(D)
    return;
end

a=D.avghet; a(a>0.5)=0.5; D.avghet=a;
D.maf=snp_het2maf(D.avghet);

  %  rs_id: {19x1 cell}
  %   allele: {19x1 cell}
  %   strand: [19x1 double]
  %  allelea: {19x1 cell}
  %  alleleb: {19x1 cell}
  %    freqa: [19x1 double]
  %    freqb: [19x1 double]
  %   neinuc: {19x1 cell}
  %   ancnuc: {19x1 cell}

[X]=snp_dbsnpinfo(D.rsid);
D.alleles=X.allele;

idx=[]; conseq=[];
  for k=1:length(D.chrpos)
      %k=15
      rspos=D.chrpos(k);
     try
      ancnuc=nt2int(D.alleles{k}(1));
      mutnuc=nt2int(D.alleles{k}(3));
      [typen]=mutationtype(seq,cdspos,strand,rspos,ancnuc,mutnuc);
     catch
         typen=-99;
     end
%  if (typen>0)
  %fprintf('%s\t%d\t%s/%s\t%d %d\n',...
  %    D.rsid{k},D.position(k),D.allele{k},D.alleleb{k},D.strand(k),strand);
%       fprintf(fid,'%s%d%d%d\t%d\t%s\t%f\t%f\t%f\t%f\n',...
%           geneid,chrid,exonstart,exonend,typen,D.rsid{k},...
%           D.afr_freq(k),D.eur_freq(k),D.chn_freq(k),D.fst3(k));
%  end
      ns=ns+(typen==1);
      nn=nn+(typen==2);
      ni=ni+(typen==0);
      %if typen>0
      %    idx=[idx,k];
      %    conseq=[conseq,typen];
      %end
      conseq=[conseq,typen];
  end
  D.conseq=conseq;
  

