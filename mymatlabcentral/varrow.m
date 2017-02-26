function h=varrow(x,ystart)

if length(x)>1
for I=1:length(x)
    h(I)=varrow(x(I),ystart);
end
else
    g=ishold(gca);
    hold on

    y=get(gca,'ylim');
    if isempty(ystart)
        h=arrow([x y(2)],[x y(1)]);
    else
        h=arrow([x ystart],[x y(1)]);
    end
    arrow(h,'EdgeColor','r')
    arrow(h,'FaceColor','r')
    
    if g==0
        hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end