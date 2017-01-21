function test_hapfreqxdist()

%%
chrid=12; startn=86000000; endn=89000000;
%chrid=12; startn=87000000; endn=87400000;
%chrid=9; startn=136300000; endn=136500000;
%chrid=2; startn=116600000; endn=118000000;

filename=sprintf('chr%d_%d_%d.mat',chrid,startn,endn);
if exist(filename,'file')
    load(filename);
else
    [hap1,pos1]=i_gethap(chrid,startn,endn,'CEU');
    [hap2,pos2]=i_gethap(chrid,startn,endn,'YRI');
    [hap3,pos3]=i_gethap(chrid,startn,endn,'ASN');
    [s1,px1]=snp_hapfreqxdist(hap1,21,1,pos1);
    [s2,px2]=snp_hapfreqxdist(hap2,21,1,pos2);
    [s3,px3]=snp_hapfreqxdist(hap3,21,1,pos3);
    save(filename);
end

%%
figure; 
subplot(2,1,1); 
plot(px1,s1,'g-');
hold on
plot(px2,s2,'k-');
plot(px3,s3,'r-');
hold off
%vline(pos1(1:20:end))
hline(0.5)

%%
close all;

idxc=10:40;
hap1c=hap1(:,idxc)+1;
hap2c=hap2(:,idxc)+1;
hap3c=hap3(1:120,idxc)+1;
%x=flipud(colormap('gray'));


figure;
subplot(2,3,4)
[h,d]=i_kmeanssort(hap1c);
snp_vhview(h); if d>0, hline(d+1,'g-'); end
xlabel(''); ylabel('')

subplot(2,3,5)
[h,d]=i_kmeanssort(hap2c);
snp_vhview(h); if d>0, hline(d+1,'g-'); end; xlabel(''); ylabel('');

subplot(2,3,6)
[h,d]=i_kmeanssort(hap3c);
snp_vhview(h); if d>0, hline(d+1,'g-'); end; xlabel(''); ylabel('');




figure;
h=subplot(2,3,1);
hapmat(hap1c)
title(snp_hapfreqxdist(hap1c))
set(h,'clim',[0 1])
%colorbar;

subplot(2,3,2)
hapmat(hap2c)
title(snp_hapfreqxdist(hap2c))

subplot(2,3,3)
hapmat(hap3c)
title(snp_hapfreqxdist(hap3c))
%xlabel(snp_hapfreqxdist(1+hap3(1:end,10:30)))
%sortfun=@sortrows;



%%
%{
figure;
subplot(3,1,1)
snp_vhview(sortrows(hap1)+1);
xlabel(''); ylabel('');
subplot(3,1,2)
snp_vhview(sortrows(hap2)+1);
xlabel(''); ylabel('');
subplot(3,1,3)
snp_vhview(sortrows(hap3(1:120,:))+1);
xlabel(''); ylabel('');
%}

end


function [hapdata,snppos]=i_gethap(chrid,startn,endn,popid)
    Data=[];
    load(sprintf('chr%d_%s_haplo',chrid,popid),'Data');
    load(sprintf('chr%d_%s_info',chrid,popid),'markinfo');
    idx=(markinfo.chrpos>startn&markinfo.chrpos<endn);
    hapdata=Data(:,idx);
    snppos=markinfo.chrpos(idx);
end