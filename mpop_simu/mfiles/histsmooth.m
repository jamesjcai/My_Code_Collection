function histsmooth(coeboot,barflag,xx)
if nargin<3
    xx='-';
end

if nargin<2
    barflag=0;
end

step=(max(coeboot)-min(coeboot))./500;
x=min(coeboot):step:max(coeboot);
y=normpdf(x,mean(coeboot),std(coeboot));

if barflag
    pdfplot(coeboot); 
    title('')
    xlabel('')
    grid off
    hold on; 
end
    plot(x,y,xx)
    vline(mean(coeboot)+1.96*std(coeboot),'b:')
    vline(mean(coeboot)-1.96*std(coeboot),'b:')
