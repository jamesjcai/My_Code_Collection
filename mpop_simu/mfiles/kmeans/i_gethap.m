function [hapdata,snppos]=i_gethap(chrid,startn,endn,popid)
    Data=[];
    load(sprintf('chr%d_%s_haplo',chrid,popid),'Data');
    load(sprintf('chr%d_%s_info',chrid,popid),'markinfo');
    idx=(markinfo.chrpos>startn&markinfo.chrpos<endn);
    hapdata=Data(:,idx);
    snppos=markinfo.chrpos(idx);
end