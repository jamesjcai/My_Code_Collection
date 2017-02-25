function histsmooth(data,barflag,xx,yy)
if nargin<4
    yy=1.5;
end

if nargin<3
    xx='-';
end

if nargin<2
    barflag=0;
end

data(isnan(data))=[];
stepn=(max(data)-min(data))/500;
x=min(data):stepn:max(data);
y=normpdf(x,mean(data),std(data));

if barflag
    pdfplot(data);
    title('')
    xlabel('')
    grid off
    hold on; 
end
plot(x,y,xx,'linewidth',yy)
%    vline(mean(data)+1.96*std(data),'b:')
%    vline(mean(data)-1.96*std(data),'b:')
