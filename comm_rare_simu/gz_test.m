%%%产生y数据
%%%第一步：产生x数据
clc
clear all
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
maf=sum(g)./n;
bar(maf)
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
% snp_writelinkage(geno,mark,'test');这里好像有错
snp_writelinkage(geno,mark,'aaa');

geno2=snp_readlinkage('aaa.ped','Delimiter',' ','UseACGT',true,'Noise',false);
if sum(sum(double(geno)-double(geno2)))~=0, error('xxx'); end
g012=snp_012geno(geno);
%第二步：gcta产生y数据
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
  
 [status2 ]=system('gcta.exe --bfile aaa --simu-qt --simu-causal-loci causal.snplist --simu-hsq 0.5 --simu-rep 1 --out aaa');
 
 else
     disp('Two SNPs are in strong LD.');
     return;
 end
else
    disp('No causal SNP(s).');
    return;
end
%第三步：eqtl分析找设置的causal SNPS
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
    %res2=fitlm(y,g012(:,k));
    %res2.Coefficients.pValue(2)
end

[pvalue,idx_best]=sort(pv);

fprintf('\n\n------------------------------\n');
fprintf('Most sig SNP is: #%d %s (P=%.2e)\n',idx_best(1),mark.rsid{idx_best(1)},pvalue(1));
fprintf('2nd most sig SNP is: #%d %s (P=%.2e)\n',idx_best(2),mark.rsid{idx_best(2)},pvalue(2));

fid=fopen('causal.snplist','r');
[D]=textscan(fid,'%s%f');
fclose(fid);
fprintf('  Causal SNP 1 (#%d) is: %s\n',idx,D{1}{1});
fprintf('  Causal SNP 2 (#%d) is: %s\n',idx2,D{1}{2});
fprintf('------------------------------\n');

%%
%%%第四步:LASSO寻找设置的casaul SNPS
% addpath('../mymatlabcentral/');
figure;
expr=y;
% evqtlplot(expr,geno(:,idx_best(1)));%%%这里好像有错，
evqtlplot(expr,geno(:,idx));
xlabel(D{1}{1})

figure;
% evqtlplot(expr,geno(:,idx_best(2)));%%%这里好像有错，
evqtlplot(expr,geno(:,idx2));
xlabel(D{1}{2})

X=g012;y=expr;
[Bb FitInfo] = lassoglm(X,y,'normal','CV',10);
minpts = find(Bb(:,FitInfo.IndexMinDeviance))
min1pts = find(Bb(:,FitInfo.Index1SE))
lassoPlot(Bb,FitInfo,'plottype','CV');
%%
figure
ypre=FitInfo.Intercept(FitInfo.Index1SE)+X*Bb(:,FitInfo.Index1SE);
 c = linspace(1,10,length(ypre));
scatter(ypre,expr,25,c,'filled')
ylabel('y predicted')
xlabel('y observed')
%%
%%%测试

[maf1,I] = sort(maf,2,'ascend'); 
g0121=g012(:,I);
nn=length(I);


B=[0.6*abs(log10(maf1(1:5)))';zeros(5,1)];
expr2=g0121(:,1:10)*B;
X=g0121;y=expr2;
[Bb2 FitInfo] = lassoglm(X,y,'normal','CV',10);
minpts = find(Bb2(:,FitInfo.IndexMinDeviance))
min1pts = find(Bb2(:,FitInfo.Index1SE))
lassoPlot(Bb2,FitInfo,'plottype','CV');
bar(Bb2)

%%
%%%%测试：GCTA产生的y于线形的y一起作用
s1=find(I==[idx])
s2=find(I==[idx2])
B=[0.6*abs(log10(maf1(1:5)))';0.6*abs(log10(maf1([s1 s2])))';zeros(5,1)];
ss=[1 2 3 4 5 s1 s2 6 7 70 71 72];%%%%根据具体的s1和s2调整最后的几个数
expr2=g0121(:,ss)*B;
X=g0121(:,ss);y=expr2;
[Bb2 FitInfo] = lassoglm(X,y,'normal','CV',10);
minpts = find(Bb2(:,FitInfo.IndexMinDeviance))
min1pts = find(Bb2(:,FitInfo.Index1SE))
lassoPlot(Bb2,FitInfo,'plottype','CV');
bar(Bb2)
%%
B=[0.6*abs(log10(maf1(1:5)))';0.6*abs(log10(maf1(75:79)))';zeros(5,1)];
%sss=[1 2 3 4 5 80 81 82 83 84 6 7 8 85 86]
sss=[1 2 3 4 5 75 76 77 78 79 6 7 8 82 83]
expr2=g0121(:,sss)*B;
X=g0121(:,sss);y=expr2;
[Bb2 FitInfo] = lassoglm(X,y,'normal','CV',10);
minpts = find(Bb2(:,FitInfo.IndexMinDeviance))
min1pts = find(Bb2(:,FitInfo.Index1SE))
lassoPlot(Bb2,FitInfo,'plottype','CV');
bar(Bb2)






