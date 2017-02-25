function [nn,ns,ni]=geneid2nsynSNP(geneid)

nn=-1; ns=-1; ni=-1;
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
    
    % get SNPs
    
popid='CEU';
%s.rs_id, s.allele, h.allelea, h.alleleb, h.freqa, h.freqb, s.ancnuc,
sql='SELECT s.rs_id, s.strand, s.position, h.allelea, h.alleleb, h.freqa, h.freqb, s.neinuc, s.ancnuc FROM snpdata.snp s inner join snpdata.hapmapii h';
sql=strcat(sql, ' on s.rs_id=h.rs_id where s.allele not like ''%%-%%''');
sql=strcat(sql,... 
sprintf(' and h.pop_id=''%s'' and s.chr_id=%d and s.position >=%d and s.position <=%d order by s.position;'...
,popid,chrid,startn,endn));

conn = database('snpdata', '', '');
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
  for k=1:length(D.position)
      %k=15
  %fprintf('%s\t%d\t%s/%s\t%d %d\n',...
  %    D.rs_id{k},D.position(k),D.allelea{k},D.alleleb{k},D.strand(k),strand);
  rspos=D.position(k);
  ancnuc=nt2int(D.allelea{k});
  mutnuc=nt2int(D.alleleb{k});
  
  [typen]=mutationtype(seq,cdspos,strand,rspos,ancnuc,mutnuc);
  ns=ns+(typen==1);
  nn=nn+(typen==2);
  ni=ni+(typen==0);
  end
  