function [geno,mark]=mpopout2genomark(x,y)
geno=snp_hap2geno(uint8(x+1));
mark.pos=round(y*500000);
mark.rsid=num2cellstr(1:length(y));