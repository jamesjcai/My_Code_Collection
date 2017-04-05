fid=fopen('aaa.phen','r');
[D]=textscan(fid,'%d%d%f');
fclose(fid);

y=D{3};
geno2=snp_readlinkage('aaa.ped','Delimiter',' ','UseACGT',true,'Noise',false);
if sum(sum(double(geno)-double(geno2)))~=0, error('xxx'); end
g012=snp_012geno(geno);

pv=ones(size(g012,2),1);
for k=1:size(g012,2)
    res=regstats(y,g012(:,k),'linear');
    pv(k)=res.fstat.pval;
end

[pvalue,idx]=min(pv);

fprintf('\n\n------------------------------\n');
fprintf('Most sig SNP is: %s (P=%.2e)\n',mark.rsid{idx},pvalue);

fid=fopen('causal.snplist','r');
[D]=textscan(fid,'%s%f');
fclose(fid);
fprintf('  Causal SNP is: %s\n',D{1}{1});
fprintf('------------------------------\n');

