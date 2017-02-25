function [bar2]=baroverlap(d,dtxt,pos)

if nargin<3
    pos=1;
end

if nargin<2
    dtxt={'a','b'};
end

[n,m]=size(d);

x=1:n;

K=1;
y1=d(:,1);
bar1=bar(x, y1, 'FaceColor', [0.50196     0.50196     0.50196], 'EdgeColor', [0.50196     0.50196     0.50196]); 
set(bar1,'BarWidth',K);
hold on;
y2=d(:,2);
bar2=bar(x, y2, 'FaceColor', 'r', 'EdgeColor', 'r');
set(bar2,'BarWidth',K*0.7); 
hold off; 

%legend(dtxt,'Interpreter','latex',pos)