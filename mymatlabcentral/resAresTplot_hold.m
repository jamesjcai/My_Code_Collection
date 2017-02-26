function [a,a2,redpoints]=resAresTplot_hold(resA,resT,ttext,isq,showcloud,allowoutlier)
if nargin<6
    allowoutlier=1;
end
if nargin<5
    showcloud=1;
end
if nargin<4
    isq=0;
end

%outidx=outlier(resA,0.005,20);
%if ~(isempty(outidx))
%resA(outidx)=[];
%resT(outidx)=[];
%end
if (showcloud)
	%plot(resA,resT,'.','markersize',1,'MarkerEdgeColor',[0.50196     0.50196     0.50196])
    colormat=[1     0.69412     0.39216];
    %colormat=[0.87059      0.4902           0];
    
    plot(resA,resT,'.','markersize',6,'MarkerEdgeColor',colormat)
    
    %plot(resA,resT,'.','markersize',1)
    %addreg(resA,resT,0,'r-')
%    addreg(resA,resT)

hold on
end



	mid=prctile(resA,[1:100]);
	[num,Xpool] = histc(resA,mid,1);

	[pooledmean]=grpstats(resT,Xpool,@median);
	[pooledmean2]=grpstats(resA,Xpool,@median);    
	uniX=unique(Xpool);
	pooledmean=pooledmean(uniX>0);
   	pooledmean2=pooledmean2(uniX>0);
	uniX=uniX(uniX>0);

	midx=mid(uniX);
%	plot(midx,pooledmean,'or','markersize',4,'MarkerFaceColor','w')
	plot(pooledmean2,pooledmean,'ko','markersize',8,...
        'MarkerFaceColor','none','LineWidth',2)	
    %plot(midx,pooledmean,'or','MarkerFaceColor','w')
    
    if nargout>2
        redpoints=[pooledmean2,pooledmean];
    end

    
    %xlim([0 midx(end-1)+1])   
 if ~allowoutlier
    xlim([0 pooledmean2(end-1)])
 end

    %if ~isq
    %ylim([0.1 .3])
    %end
	axis square
	hold off


[a,p1]=corr(midx(1:end-1)',pooledmean(1:end-1),'type','spearman');
[a2,p2]=corr(resA,resT,'type','spearman');

%[a,p1]=corr(midx(1:end-1)',pooledmean(1:end-1));    
%if isq
%    fprintf('Qneu\t%f (%.4f), %f (%.4f) %s\n',a,p1,a2,p2,ttext)       
%else


%fprintf('Pneu\t%f (%.4f), %f (%.4f) %s\n',a,p1,a2,p2,ttext)

 %{   
   if p2<0.001
        fprintf('%s\t%.4f (p=%.4f)\n',ttext,a2,p2);
        xxa=sprintf('r=%.4f^{**}\n',a2);        
    else
        fprintf('%s\t%.4f (%.4f)',ttext,a2,p2);
        xxa=sprintf('r=%.4f (p=%.4f)',a2,p2);                
    end
   %} 
    
%    legend(xxa)
%    legend('boxoff')

if p2<0.05
    tag2='^*';
else
    tag2='^{ns}';
end
if p1<0.05
    tag1='^*';
else
    tag1='^{ns}';
end

xlabel(sprintf('%s\n\nr=%.3f%s, r_{pooled}=%.3f%s',ttext,a2,tag2,a,tag1))
%xlabel(sprintf('%s',ttext))
if isq
	ylabel('Qneu')    
else
	ylabel('Pneu')
end

%h=gcf;set(h,'FontName','Arial')
