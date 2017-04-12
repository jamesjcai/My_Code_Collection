[status1]=system('plink --file aaa --make-bed --out aaa --silent');


idxv=find(maf>0.2);   % common SNP list
idxvr=find(maf<0.05); % rare SNP list

if isempty(idxv) || length(idxv)<1 || length(idxvr)<1
    disp('No causal SNP(s).');
    return;
end
  
 idx1=idxv(1);
 idx2=idxv(end);
 
 dstat=snp_ldpair(snp_pickmarker(geno,[],[idx1,idx2]));
 dstat.dprime 
 if dstat.dprime(1,2)>0.85, disp('Two SNPs are in strong LD.'); return; end
 
  
 
 fid=fopen('common.causal.snplist','w');
 fprintf(fid,'%s\t-0.025\n',mark.rsid{idx1});
 fprintf(fid,'%s\t0.021\n',mark.rsid{idx2});
 % fprintf(fid,'%s\t3.0\n',mark.rsid{idx3});
 fclose(fid); 
 
 [status2]=system('gcta.exe --bfile aaa --simu-qt --simu-causal-loci common.causal.snplist --simu-hsq 0.5 --simu-rep 1 --out aaa');
 % system('gcta --bfile aaa --simu-cc 500 500 --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-k 0.1 --simu-rep 3 --out aaatest');

fid=fopen('aaa.phen','r');
[D]=textscan(fid,'%d%d%f');
fclose(fid);
expr_common=D{3};




 idx3=idxvr(randperm(length(idxvr)));
 idx3=idx3(1:round(length(idx3)*0.1));
 fprintf('\n\nAdding effect of %d rare SNPs...\n',length(idx3));
 
 geno_rare=snp_pickmarker(geno,[],idx3);
 g012_rare=snp_012geno(geno_rare); 
 % g012_rare2=g012(:,idx3);
 % if sum(sum(double(g012_rare2)-double(g012_rare)))~=0, error('xxx1'); end

 % B=randn(length(idx3),1).*(1000*maf(idx3).^2)';
 B = betapdf(maf(idx3),1,25)';
 
 expr = expr_common;% + g012_rare*B;
 % expr = g012_rare*B;
 
 
 
  