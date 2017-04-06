[status1]=system('plink --file aaa --make-bed --out aaa --silent');

idxv=find(maf>0.3);

if ~isempty(idxv) && length(idxv)>1
 % idx=idx(randperm(length(idx)));
 idx=idxv(1);
 idx2=idxv(end);
 
 fid=fopen('causal.snplist','w');
 fprintf(fid,'%s\t-0.025\n',mark.rsid{idx});
 fprintf(fid,'%s\t0.021\n',mark.rsid{idx2});
 fclose(fid);
 
 dstat=snp_ldpair(snp_pickmarker(geno,[],[idx,idx2]));
 dstat.dprime
 
 if dstat.dprime(1,2)<0.8    
  
 [status2]=system('gcta.exe --bfile aaa --simu-qt --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-rep 1 --out aaa');
 % system('gcta --bfile aaa --simu-cc 500 500 --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-k 0.1 --simu-rep 3 --out aaatest');
 else
     disp('Two SNPs are in strong LD.');
     return;
 end
else
    disp('No causal SNP(s).');
    return;
end

