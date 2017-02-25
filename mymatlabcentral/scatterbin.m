function [a,a2,redpoints]=scatterbin(x,y,showcloud,colorset,allowoutlier)

%SEE also: pooledplot (mean)
if nargin<5, allowoutlier=true; end
if nargin<4, colorset=1; end
if nargin<3, showcloud=true; end

%{
x=x(~isnan(y));
y=y(~isnan(y));
y=y(~isnan(x));
x=x(~isnan(x));
%}
idx=~(isnan(x)|isnan(y));
x=x(idx);
y=y(idx);

if (showcloud)
    %plot(x,y,'.','markersize',1,'MarkerEdgeColor',[0.50196 0.50196 0.50196])    
    %plot(x,y,'.','markersize',1,'MarkerEdgeColor','w')
    if colorset==2
        plot(x,y,'+','markersize',3,'MarkerEdgeColor',[1 0.6 0.78431]);
        %i_addreg(x,y,1,'r-',1);        
    elseif colorset==1
        plot(x,y,'x','markersize',3,'MarkerEdgeColor',[0.70196     0.78039           1]);
        %i_addreg(x,y,1,'b-',1);        
    end
    hold on
end

	mid=prctile(x,1:100);
	[~,Xpool] = histc(x,mid,1);

	[pooledmean]=grpstats(y,Xpool,@median);
	[pooledmean2]=grpstats(x,Xpool,@median);    
	uniX=unique(Xpool);
	pooledmean=pooledmean(uniX>0);
   	pooledmean2=pooledmean2(uniX>0);
	uniX=uniX(uniX>0);

	midx=mid(uniX);
    %	plot(midx,pooledmean,'or','markersize',4,'MarkerFaceColor','w')
    if colorset==1    
        plot(pooledmean2,pooledmean,'or','markersize',5,'MarkerFaceColor','w')
    elseif colorset==2
        plot(pooledmean2,pooledmean,'sb','markersize',5,'MarkerFaceColor','w')
    end
    
    if colorset==1
        i_addreg(pooledmean2,pooledmean,1,'r-',1);
    elseif colorset==2
        i_addreg(pooledmean2,pooledmean,1,'b-',1);
    end
    
    
    %   plot(midx,pooledmean,'or','MarkerFaceColor','w')    
    if nargout>2
        redpoints=[pooledmean2,pooledmean];
    end
    
    %xlim([0 midx(end-1)+1])   
 if ~allowoutlier
    xlim([0 pooledmean2(end-1)])
 end
%	axis square
	hold off


[a,p1]=corr(midx(1:end-1)',pooledmean(1:end-1),'type','spearman');
[a2,p2]=corr(x,y,'type','spearman');

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

if p2<0.05, tag2='^*'; else tag2='^{ns}'; end

if p1<0.05, tag1='^*'; else tag1='^{ns}'; end

% title(sprintf('Spearman rho=%.3f%s, rho_{pooled}=%.3f%s',a2,tag2,a,tag1))
title(sprintf('r=%.2f^{ p=%.1e}, r_{pooled}=%.2f^{ p=%.1e}',a2,p2,a,p1))
%h=gcf;set(h,'FontName','Arial')



function i_addreg(Dn,Pneu,showqudratic,linestyle,lineweight)
%i_addreg(Dn,Pneu,showqudratic,linestyle)

if nargin<5
    lineweight=2;
end

if nargin<4
    linestyle='r-';
end

if nargin<3
    showqudratic=false;
end

Dn=Dn(~isnan(Pneu));
Pneu=Pneu(~isnan(Pneu));
Pneu=Pneu(~isnan(Dn));
Dn=Dn(~isnan(Dn));


%plot(Dn, Pneu,'.')
[orixy]=axis;

x=[min(Dn),max(Dn)];
[a]=regress1d(Dn, Pneu);
y=a(1)*x+a(2);
hold on;

plot(x,y,linestyle,'LineWidth',lineweight)
xlim([orixy(1) orixy(2)])
ylim([orixy(3) orixy(4)])

%{
try
[a]=L1LinearRegression(Dn, Pneu);
catch
[a]=L1LinearRegression(Dn', Pneu');
end
y=a(1)*x+a(2);
plot(x,y,'b-')
%}

%[r,p]=corr(Dn,Pneu,'type','spearman');
%fprintf('Spearman rho=%r (P=%f)\n',r,p);

if showqudratic
    % Calculate fit parameters
    p=polyfit(Dn,Pneu,2);
    x2=linspace(min(Dn),max(Dn));
    y2=polyval(p,x2);
    %plot(x2,y2,'--','color',[0.25098     0.50196     0.50196]);
    %plot(x2,y2,'g--','LineWidth',2);
    %plot(x2,y2,linestyle,'LineWidth',lineweight,'linesmoothing','on');
    plot(x2,y2,linestyle,'LineWidth',lineweight);
end
hold off;

