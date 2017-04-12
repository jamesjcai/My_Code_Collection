


%geno2=snp_readlinkage('aaa.ped','Delimiter',' ','UseACGT',true,'Noise',false);
%if sum(sum(double(geno)-double(geno2)))~=0, error('xxx'); end

pv=ones(size(g012,2),1);
for k=1:size(g012,2)
    res=regstats(expr,g012(:,k),'linear');
    pv(k)=res.fstat.pval;
    %res2=fitlm(y,g012(:,k));
    %res2.Coefficients.pValue(2)
end

[pvalue,idx_best]=sort(pv);

fprintf('\n\n------------------------------\n');
fprintf('Most sig SNP is: #%d %s (P=%.2e)\n',idx_best(1),mark.rsid{idx_best(1)},pvalue(1));
fprintf('2nd most sig SNP is: #%d %s (P=%.2e)\n',idx_best(2),mark.rsid{idx_best(2)},pvalue(2));
fprintf('3rd most sig SNP is: #%d %s (P=%.2e)\n',idx_best(3),mark.rsid{idx_best(3)},pvalue(3));

fid=fopen('common.causal.snplist','r');
[D]=textscan(fid,'%s%f');
fclose(fid);
fprintf('  Causal SNP 1 (#%d) is: %s\n',idx1,D{1}{1});
fprintf('  Causal SNP 2 (#%d) is: %s\n',idx2,D{1}{2});
if length(D{1})>2
    fprintf('  Causal SNP 3 (#%d) is: %s\n',idx3,D{1}{3});
end
fprintf('------------------------------\n');

