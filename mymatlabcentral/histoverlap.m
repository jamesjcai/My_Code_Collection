function histoverlap(data1,data2,binnum,legtxt,logtag)
% URL: http://desk.stinkpot.org:8080/tricks/index.php/category/matlab/

if nargin<1, data1=randn(1,1000); end
if nargin<2, data2=randn(1,1000)+2; end
if nargin<3, binnum=30; end
if nargin<4, legtxt={'A','B'}; end
if nargin<5, logtag=false; end

if ~logtag
    hist(data1,binnum)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
    hold on
    hist(data2,binnum)
    h = findobj(gca,'Type','patch');
    set(h,'facealpha',0.75);
    ylabel('counts')
    %xlabel('gene-tree-length/species-tree-length')
    legend(legtxt,2)
    %title('p = 0.5');

else    
    [y,x]=hist(data1,binnum);
    h=bar(x,y);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
    hold on
    [y,x]=hist(data2,binnum);
    h=bar(x,y);
    h = findobj(gca,'Type','patch');    
    set(h,'facealpha',0.75);
    ylabel('counts')
    set(gca,'yscale','log')
        legend(legtxt)

end
%set(h,'BaseValue',0.1)


%legend([repmat('penalty = ',length(vec),1) num2str(vec)])