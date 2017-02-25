function [D]=perlegen_at_region(chrid,startn,endn)


	%mmFilename=sprintf('C:/biodata/hgenome/Homo_sapiens.NCBI36.40.dna.chromosome.%d.mm',chrid);
	%chrmm = memmapfile(mmFilename, 'format', 'uint8');
        %seq=chrmm.Data(startn:endn)';

%get SNPs

%popid='CEU';
%s.rs_id, s.allele, h.allelea, h.alleleb, h.freqa, h.freqb, s.ancnuc,

sql='SELECT p.rsid, p.chromosome, p.position, p.alleles, p.afr_freq, p.eur_freq, p.chn_freq, p.fst_ea, p.fst_ac, p.fst_ec, p.fst3';
sql=strcat(sql, ' FROM Perlegenchrall p');
sql=strcat(sql,...
sprintf(' WHERE (((p.chromosome)=%d) AND ((p.position)>=%d) AND ((p.position)<%d));',...
        chrid,startn,endn));

conn = database('perlegen', '', '');
curs = exec(conn, sql);
setdbprefs('DataReturnFormat','structure');
curs = fetch(curs);
D=curs.Data;
close(curs);
close(conn);

if ~isstruct(D)
    D=[];
end

%    D.chrid=chrid;
%    D.startn=startn;
%    D.endn=endn;
