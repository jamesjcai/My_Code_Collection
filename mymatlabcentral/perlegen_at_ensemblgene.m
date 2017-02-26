function [nn,ns,ni,D]=perlegen_at_ensemblgene(geneid)

nn=-1; ns=-1; ni=-1; D=[];
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

%popid='CEU';
%s.rs_id, s.allele, h.allelea, h.alleleb, h.freqa, h.freqb, s.ancnuc,

sql='SELECT p.rsid, p.chromosome, p.position, p.alleles, p.afr_freq, p.eur_freq, p.chn_freq, p.fst_ea, p.fst_ac, p.fst_ec, p.fst3';
sql=strcat(sql, ' FROM Perlegenchrall p');
sql=strcat(sql,...
sprintf(' WHERE (((p.chromosome)=%d) AND ((p.position)>=%d) AND ((p.position)<%d));',...
        chrid,startn,endn));


conn = database('perlegen', '', '');
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
      ancnuc=nt2int(D.alleles{k}(1));
      mutnuc=nt2int(D.alleles{k}(3));
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

if (nargout>3 && ~isempty(idx))
	D.rsid=D.rsid(idx);
	D.alleles=D.alleles(idx);
	D.position=D.position(idx);
	D.afr_freq=D.afr_freq(idx);
	D.eur_freq=D.eur_freq(idx);
	D.chn_freq=D.chn_freq(idx);
    
	D.fst_ea=D.fst_ea(idx);
    D.fst_ac=D.fst_ac(idx);
    D.fst_ec=D.fst_ec(idx);
    D.fst3=D.fst3(idx);
    
	D.conseq=conseq;
	D.chromosome=chrid;
	D.ancnuc=chrnuc_chimp(chrid,D.position);
	for (k=1:length(D.ancnuc)),
	      ancx=D.ancnuc(k);
	      mutx1=nt2int(D.alleles{k}(1));
	      mutx2=nt2int(D.alleles{k}(3));
	      if (ancx==mutx1)
		     D.daf(k)=1-D.chn_freq(k);
	      elseif (ancx==mutx2)
		     D.daf(k)=D.chn_freq(k);
          else
		    D.daf(k)=-1;
	      end
	end
else
	D=[];
end



