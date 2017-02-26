function genometrackx(chrid,regionl,regionr)

if nargin==0
    chrid=17;
    regionl=43e6;
    regionr=46e6;
end
load namedgenes_hs_ncbi_b36;

%xls_chrid	xls_strand	xls_genestart	xls_geneend	xls_ensembid
%xls_genename
idx=(xls_chrid==chrid&xls_genestart>=regionl&xls_geneend<=regionr);

xgenename=xls_genename(idx);
%xstrand=xls_strand(idx);
xgenestart=xls_genestart(idx);
xgeneend=xls_geneend(idx);
if ~isempty(xgenestart)
    i_geneplot([xgenestart xgeneend],xgenename,[regionl,regionr])
else
    plot([regionl,regionr],[0 0],'-k');
end
xlim([regionl,regionr])

function i_geneplot(params,gene,regions)

    %plot(regions./1000000,[0 0],'-k');
    %params=params./1000000;
    plot(regions,[0 0],'-k');
    nb = size(params,1);  %  number of boxplots
    hw = 0.15;   %  halfwidth of boxes

    %plot([min(min(params)) max(max(params))],[0 0],'k-')
    c=-1;  
    hold on
    for ii = 1:nb
       temp1 = params(ii,1);
       temp2 = params(ii,2);
       xx = [temp1 temp1 temp2 temp2 temp1];
    %   yy = [ii-hw ii+hw ii+hw ii-hw ii-hw];
       yy = [-hw hw hw -hw -hw];
       %plot(xx,yy,'-')
       %fill(xx,yy,'r');
       patch(xx,yy,'w');
       %text(temp1,hw,gene{ii},'fontsize',7,'Rotation',90);
       c=c*-1;
       if c>0
        text(temp1+round(0.5*(temp2-temp1)),hw+(hw*0.3),gene{ii},'fontsize',7,'Rotation',90);
       else
        text(temp1+round(0.5*(temp2-temp1)),-hw-(hw*0.3),gene{ii},'fontsize',7,'Rotation',-90);
       end
    end
    hold off

    %  make some extra space
    %axlim = axis;
    %axlim(1) = axlim(1)-1;
    %axlim(2) = axlim(2)+1;
    %axlim(3) = axlim(3)-2;
    %axlim(4) = axlim(4)+2;
    %axis(axlim)
    
    ylim([-2 2])
    %{
    [a]=get(gca,'position');
    a(2)=0.5;
    a(4)=0.3;
    set(gca,'position',a);
    %}
    set(gca,'yticklabel',[],'ytick',[]);
    %set(gca,'xcolor',get(clf,'color'),'xtick',[]);


