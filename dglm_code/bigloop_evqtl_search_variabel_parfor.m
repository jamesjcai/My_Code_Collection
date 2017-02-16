function bigloop_evqtl_search_variabel_parfor(z1,z2)
if nargin<2
    error('need two input values: e.g., program 1000 2000\n');
end
if isdeployed
    z1=str2double(z1);
    z2=str2double(z2);
end

warning('off','all')
%{
if isunix
load ('/media/cailab/Seagate3TB1/home/cailab/NetDrive3TB/gwang/wg/1000_genome_evQTL/gene_exp/PeerRPKM_406_no_outliers_protein_coding');
load ('/media/cailab/Seagate3TB1/home/cailab/NetDrive3TB/gwang/wg/1000_genome_evQTL/gene_location/human_GRCh37_gene_location');
else
%load ('\\165.91.29.77\NetDrive3TB\gwang\wg\1000_genome_evQTL\gene_exp\PeerRPKM_406_no_outliers_protein_coding');
%load ('\\165.91.29.77\NetDrive3TB\gwang\wg\1000_genome_evQTL\gene_location\human_GRCh37_gene_location');
end

%%
idx=~ismember(groupid,'AFR');
X6=X5(:,idx);
idv6=idv5(idx);

%%
if isunix
load ('/media/cailab/Seagate3TB1/home/cailab/NetDrive3TB/gwang/wg/1000_genome_evQTL/snp_geno/snp_1000_all_chr_005');
else
%load ('\\165.91.29.77\NetDrive3TB\gwang\wg\1000_genome_evQTL\snp_geno\snp_1000_all_chr_005')
end

snp_geno6=snp_geno(:,idx);

%}
load bigloop_evqtl_search_data.mat

%load logfile.txt
%currentid=max(logfile)+1;


%%
if z2>size(X6,1), z2=size(X6,1); end
parpool
for k=z1:z2 % currentid:size(X6,1)
    %fid2 = fopen('logfile.txt','a');

    tic
    %fprintf(fid2,'Start gene %d\t%s......\n',k,geneid5{k});
    fprintf('gene %d\t%s......\n',k,geneid5{k});
    parfor l=1:size(snp_geno6,1)
        fid = fopen(sprintf('Result_%d_%d.txt',z1,z2),'a');       
        x=snp_geno6(l,:);
        if length(unique(x))==3
            y=X6(k,:);    
            [p,effsize]=run_variABEL(zscore(y'),double(x'));
            if p<=1e-6
                fprintf(fid,'%d\t%d\t%e\t%f\n',k,l,p,effsize);
                fprintf('%d\t%d\t%e\t%f\n',k,l,p,effsize);
            end
	       end
        fclose(fid); 
    end
    %fprintf(fid2,'End gene %d\t%s......\n',k,geneid5{k});
    %fprintf(fid2,'%d\n',k);
    toc

	%fclose(fid2);
end
