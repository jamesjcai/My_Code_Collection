function [pslope1]=pooledplot(resA,resT,drawcloud)

% SEE ALSO: SCATTERBIN (median)
pslope1=1;
if (nargin<3), drawcloud=1; end
idx=~(isnan(resA)|isnan(resT));
resA=resA(idx);
resT=resT(idx);
%scatter(resA,resT,'.');
%return;


    if drawcloud
      plot(resA,resT,'.','markersize',1,'color',[0.70196     0.78039           1]);
      hold on
    end
    
	mid=prctile(resA,[1:100]);
	[num,Xpool] = histc(resA,mid,1);
	[pooledmean]=grpstats(resT,Xpool,@mean);
	uniX=unique(Xpool);
	pooledmean=pooledmean(uniX>0);
	uniX=uniX(uniX>0);
	midx=mid(uniX);
    
        midx(1)=[];
        pooledmean(1)=[];  
        midx(end)=[];
        pooledmean(end)=[];  
        
	plot(midx,pooledmean,'or','markersize',4,'MarkerFaceColor','w')
    if (1)
        %midx(end)=[];
        %pooledmean(end)=[];  
%        ylim([min(pooledmean) max(pooledmean)])
%        xlim([min(midx) max(midx)])
      
      pslope=polyfit(midx',pooledmean,1);        
      plot([0 max(midx)],[pslope(2),max(midx)*pslope(1)+pslope(2)],'g');
      pslope1=pslope(1);
      
       % ylim([0 max(pooledmean)])
       % xlim([0 max(midx)])      
    else
        midx(end)=[];
        pooledmean(end)=[];
        z=max(max(midx),max(pooledmean));
        xlim([0 z])
        ylim([0 z])
        line
        axis square
    end
    
	hold off