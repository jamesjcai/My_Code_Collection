clc
clearvars
close all
system('ms.exe 2000 1 -t 16.4 -r 100.0 2501 -G 6.93 -eG 0.2 0.0 -eN 0.3 0.5 >ms_output.txt');

OUT=readmsoutput('ms_output.txt');
gametes=(OUT.gametes{1});
[n,p]=size(gametes);
gametes=gametes(randperm(n),:);
g=zeros(n/2,p);
for k=1:2:n
 g((k+1)/2,:)=sum(gametes([k,k+1],:));
end
mafv=sum(g)./n;
mafv(mafv>0.5)=1-mafv(mafv>0.5);

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

%geno2=snp_readlinkage('aaa.ped','Delimiter',' ','UseACGT',true,'Noise',false);
%if sum(sum(double(geno)-double(geno2)))~=0, error('xxx'); end
g012=snp_012geno(geno);

% genotype simulation is done.
%%

status1=1; status2=1;
i_common_run_gcta

%%
if status1==0 && status2==0
    i_common_detect_eqtl
    i_common_lasso
end
return;

%%
save skatinfile g012 expr expr_residual
% saveR('skatinfile.r', 'g012', 'expr')
[status3]=system('E:/R/R-3.3.1/bin/x64/Rscript.exe run_SKAT.r');
if status3==0
load('skatoutfile.mat','pskat1','pskat2')
pskat1
pskat2
end





