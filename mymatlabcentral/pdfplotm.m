function pdfplotm(datalist)



linest={'b','r','g','m','c','k','g','y','r','m','c','b','k','g','y','r'};
mlinest={'.','o','x','+','*','s','d','v','^','<','>','p','h'};

h2=[];
hold on
for k=1:length(datalist)
    [f,xi] = ksdensity(datalist{k});

    %h=pdfplot(datalist{k});
    %set(h,'FaceColor','none','EdgeColor',[0.75294     0.75294     0.75294]);    
            %[0.50196     0.50196     0.50196])
    h2=[h2,plot(xi,f,sprintf('%s%s',linest{k},mlinest{k}))];
    
end
legend(h2,mat2cellstr(1:length(datalist)));
hold off
