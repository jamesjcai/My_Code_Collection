function evqtlplot(y,g,colorid)
if nargin<3
    colorid=1;
end
    c1=[0.70196 0.78039 1]; 
    c2=[1 0.6 0.78431]; 
    c3=[.8 .8 .8];
    c2=defaultcolor(2);
 switch colorid
     case 1
         c=c1;
     case 2
         c=c2;
     case 3
         c=c3;
     otherwise
         c='b';
 end
          %defaultc=get(gca,'ColorOrder');
          defaultc=get(groot,'DefaultAxesColorOrder');
          c=defaultc(colorid,:);
        %plot(1+g+0.1*(rand(size(g))-0.5),y,'o','linesmooth','off','color',c)
        %plot(1+g+0.1*(rand(size(g))-0.5),y,'o','color',c)
        
        %hold on
        h2=boxplot(y,g,'colors','k');
        set(h2(6,:),'color','k','linewidth',2.5);
        
        hold on
        %pp1=plot(1+g+0.1*(rand(size(g))-0.5),y,'o','color',c);
        pp1=scatter(1+g+0.1*(rand(size(g))-0.5),y);
        %pp1.MarkerFaceColor=c;
        pp1.MarkerEdgeColor=c;
        %pp1.MarkerFaceAlpha=0.3;
        
        %set(h2,'linesmooth','off','linewidth',1.5);

        box on
        hold off
        
        % http://stackoverflow.com/questions/21999451/how-to-get-the-values-of-the-outliers-and-their-coordinates-from-a-box-plot
