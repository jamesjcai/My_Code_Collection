function [nn,ns,ni,D]=watsonsnp_at_ensemblexon(geneid,transid,exonstart,exonend,logfile)

if nargin<5
    logfile=0;
end

nn=-1; ns=-1; ni=-1;
%{
transid='ENST00000360131';
geneid='ENSG00000084693';
genename='AGBL5';
exonstart=27129774;
exonend=27146635;
%}

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

%popid='CEU';
%s.rs_id, s.allele, h.allelea, h.alleleb, h.freqa, h.freqb, s.ancnuc,

sql='SELECT p.rsid, p.chromosome, p.position, p.alleleA, p.alleleB, p.Coverage, p.hethom';
sql=strcat(sql, ' FROM [Watson-454-snp-v01] p');
sql=strcat(sql,...
sprintf(' WHERE (((p.chromosome)=''chr%d'') AND ((p.position)>=%d) AND ((p.position)<%d));',...
        chrid,startn,endn));

conn = database('watsonsnp', '', '');
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

  %  rs_id: {19x1 cell}
  %   allele: {19x1 cell}
  %   strand: [19x1 double]
  %  allelea: {19x1 cell}
  %  alleleb: {19x1 cell}
  %    freqa: [19x1 double]
  %    freqb: [19x1 double]
  %   neinuc: {19x1 cell}
  %   ancnuc: {19x1 cell}


  ns=0;
  nn=0;
  ni=0;

%  fid=fopen('c:/perlegen_at_ensemblexon_log.txt','w');

  idx=[]; conseq=[];
  for k=1:length(D.position)
      %k=15
      rspos=D.position(k);
      ancnuc=nt2int(D.alleleA{k});
      mutnuc=nt2int(D.alleleB{k});
      [typen]=mutationtype(seq,cdspos,strand,rspos,ancnuc,mutnuc);
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
      if typen>0
          idx=[idx,k];
          conseq=[conseq,typen];
      end
  end
%fclose(fid);

if (nargout>3 & ~isempty(idx))
	D.rsid=D.rsid(idx);
	D.alleleA=D.alleleA(idx);
    D.alleleB=D.alleleB(idx);    
	D.position=D.position(idx);	
    D.coverage=D.Coverage(idx);
	D.hethom=D.hethom(idx);
    
	D.conseq=conseq;
	D.chromosome=chrid;
%	D.ancnuc=chrnuc_chimp(chrid,D.position);

else
	D=[];
end
