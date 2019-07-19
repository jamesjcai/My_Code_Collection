if ispc
    addpath('C:\Users\jcai\Documents\GitHub\My_Code_Collection\MrGSEA_rank\GSEA_package');
    
else
    addpath('/media/cailab/DISK4T/GitHub/My_Code_Collection/MrGSEA_rank/GSEA_package');
    addpath('/mnt/DISK4T/GitHub/My_Code_Collection/MrGSEA_rank/GSEA_package');
    addpath('/mnt/RD_TERRA3/DATA/DISK4T/ref_gene_sets/');
end

opts = default_GSEA_opts();
opts.show = true;     % if plot results
opts.save = true;      % if save results
opts.perm_nb = 1000;    %number of permutations

T=readtable('aaa.txt');
genelist=string(T.glist01);
generank=T.dd;

if ispc
    load \\cvm-research-dr.cvm.tamu.edu\CaiLab\Cai-Terra3\DATA\DISK4T\ref_gene_sets\msigdb_v62\msigdb_v62_c5.mat
else
    load /mnt/RD_TERRA3/DATA/DISK4T/ref_gene_sets/msigdb_v62/msigdb_v62_c5.mat
end

[res_pos,res_neg,res_descr,p_gene] = MrGSEA_rank(generank,genelist,...
                                     GeneSet,GeneSetName,'test_rank2',opts);

                                 
% load test_rank2                                 
Tp=cell2table(res_pos);
Tp.Properties.VariableNames=res_descr(1:end-1);
Tp = sortrows(Tp,'NES_qval','ascend');
Tp.GSet=GeneSetName(Tp.GSNAME);
Tn=cell2table(res_neg);
Tn.Properties.VariableNames=res_descr(1:end-1);
Tn = sortrows(Tn,'NES_qval','ascend');
Tn.GSet=GeneSetName(Tn.GSNAME);

