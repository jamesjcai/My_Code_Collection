geno=snp_readlinkage('mypedfile_tab.ped','Delimiter','\t','MissingGenotype',...
                     '0','UseACGT',true,'Noise',true);
X=snp_012geno(geno);
clear geno;
