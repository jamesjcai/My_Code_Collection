
load('\\165.91.29.187\disk4t\1000GenomeGenotype\phase_3\panel2504.mat', 'super_pop')
load('\\165.91.29.187\disk4t\1000GenomeGenotype\phase_3\mat\geno_phase3_v5a_20130502_chr22.mat')
load('\\165.91.29.187\disk4t\1000GenomeGenotype\phase_3\mat\marklite_phase3_v5a_20130502_chr22.mat')

%%
geno=geno(ismember(super_pop,'EUR'),:);
idx=pos>20000000 & pos<20300000;
[geno]=snp_pickmarker(geno,[],idx);
mark.rsid=rsid(idx);
mark.pos=pos(idx);
mark.chrid=22*ones(size(mark.pos));

snp_writelinkage(geno,mark,'aaa');
maf=snp_maf(geno);


status1=1; status2=1;
i_common_run_gcta

if status1==0 && status2==0
    i_common_detect_eqtl    
    i_common_lasso
end

