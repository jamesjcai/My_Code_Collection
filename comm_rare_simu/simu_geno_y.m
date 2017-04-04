system('ms 2000 1 -t 6.4 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5 >ms_output_example.txt');

OUT=readmsoutput('ms_output_example.txt');
gametes=(OUT.gametes{1});
[n,p]=size(gametes);
gametes=gametes(randperm(n),:);
g=zeros(n/2,p);
for k=1:2:n
 g((k+1)/2,:)=sum(gametes([k,k+1],:));
end
maf=sum(g)./n;
maf(maf>0.5)=1-maf(maf>0.5);

mark.chrid=ones(p,1);
mark.pos=round(OUT.positions{1}*10000);
mark.rsid=cell(p,1);

geno=zeros(n/2,p*2);
AC=[1 1; 1 2; 2 2];
for k=1:p
 geno(:,[k*2-1 k*2])=AC(1+g(:,k),:);
 mark.rsid{k}=sprintf('mrk%d',k);
end


snp_writelinkage(geno,mark,'aaa');

system('plink --file aaa --make-bed --out aaa');

idx=find(maf>0.3);
if ~isempty(idx)
 idx=idx(1);
 fid=fopen('causal.snplist','w');
 fprintf(fid,'mrk%d\t0.125\n',idx);
 fclose(fid);
 
 system('gcta --bfile aaa --simu-qt --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-rep 3 --keep aaatest.indi.list --out aaatest');
 % system('gcta --bfile aaa --simu-cc 500 500 --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-k 0.1 --simu-rep 3 --out aaatest');
end
