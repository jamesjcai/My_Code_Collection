function cdfplotm(datalist)


linest={'b','r','g','m','c','k','g','y','r','m','c','b','k','g','y','r'};
%mlinest={'none','.','o','x','+','*','s','d','v','^','<','>','p','h'};
%mlinest={'none','none','o','x','+','*','s','d','v','^','<','>','p','h'};
mlinest={'none','none','none','none','none','none','none','none','none','none','none','none'};

hold on
for k=1:length(datalist)
    [h,stats]=cdfplot(datalist{k});
    set(h,'color',linest{k},'marker',mlinest{k},'LineStyle','-')
end
legend(mat2cellstr(1:length(datalist)),4);
hold off