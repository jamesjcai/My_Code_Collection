[status1]=system('plink --file aaa --make-bed --out aaa --silent');

idxvc=find(mafv>0.2);  % common SNP list
idxvr=find(mafv<0.05); % rare SNP list

if isempty(idxvc) || length(idxvc)<1 || length(idxvr)<1
    disp('No causal SNP(s).');
    return;
end
  
 idx1=idxvc(1);
 idx2=idxvc(end);
 
 dstat=snp_ldpair(snp_pickmarker(geno,[],[idx1,idx2]));
 dstat.dprime 
 if dstat.dprime(1,2)>0.85, disp('Two SNPs are in strong LD.'); return; end

 fid=fopen('common.causal.snplist','w');
 fprintf(fid,'%s\t-0.25\n',mark.rsid{idx1});
 fprintf(fid,'%s\t0.21\n',mark.rsid{idx2});
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
 fprintf('Adding effect of %d rare SNPs...\n',length(idx3));
 
 geno_rare=snp_pickmarker(geno,[],idx3);
 g012_rare=snp_012geno(geno_rare);
 % g012_rare2=g012(:,idx3);
 % if sum(sum(double(g012_rare2)-double(g012_rare)))~=0, error('xxx1'); end

 % B=randn(length(idx3),1).*(1000*mafv(idx3).^2)';
 B=betapdf(mafv(idx3),1,25)'./10;
 %B=0.6*abs(log10(mafv(idx3)))';
 i=randn(length(idx3),1)<0;
 B(i)=-B(i);

 % expr = expr_common;% + g012_rare*B;
 % expr = g012_rare*B;
 expr = expr_common + g012_rare*B  + randn(length(expr_common),1);
 
 %% Display MAF and simulated Beta values. 
 betavector=zeros(length(mafv),1);
 betavector([idx1 idx2])=[-0.25 0.21];
 betavector(idx3)=B;
 
 fb=figure;
 subplot(6,1,1)
     bar(mafv)
     hold on
     line([idx1 idx1],ylim,'color','r')
     line([idx2 idx2],ylim,'color','m')
     for k=1:length(idx3)
        line([idx3(k) idx3(k)],ylim,'color','g')
     end
     line(xlim,[0.05 0.05],'color','r');
     ylim([0 0.5]);
     title('MAF and causal mutations (r,m,g)')
     ylabel('MAF')
  subplot(6,1,2)
     bar(betavector)
     title('Beta values simulated')
     ylabel('Beta')