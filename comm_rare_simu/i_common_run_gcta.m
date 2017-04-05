[status1]=system('plink --file aaa --make-bed --out aaa');

idxv=find(maf>0.3);

if ~isempty(idxv)
 % idx=idx(randperm(length(idx)));
 idx=idxv(1);
 
 fid=fopen('causal.snplist','w');
 fprintf(fid,'%s\t0.125\n',mark.rsid{idx});
 fclose(fid);
 
 [status2]=system('gcta.exe --bfile aaa --simu-qt --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-rep 1 --out aaa');
 % system('gcta --bfile aaa --simu-cc 500 500 --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-k 0.1 --simu-rep 3 --out aaatest');
else
    disp('No causal SNP.');
    return;
end