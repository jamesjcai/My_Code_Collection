function i_lineto(u,v,methodid)
if nargin<3, methodid=1; end
if nargin<2
    v=u;
    u=[0;0];
end

switch methodid
    case 1
        hold on
        line([u(1), v(1)],[u(2), v(2)],'Color','k');
        plot(u(1),u(2),'.k','MarkerSize',20);
        plot(v(1),v(2),'.k','MarkerSize',10);
    case 2
        if u(1)==0 && u(2)==0
            quiver(0,0,v(1),v(2));
        else
            error('Only for vector starting from [0 0]');
        end
end

end
