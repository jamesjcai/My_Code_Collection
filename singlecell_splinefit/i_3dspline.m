function [d,xyz1]=i_3dspline(x,y,z)

% see usage: s_test_example.m
% scatter3(x,y,z,'o','MarkerFaceColor',[0 .75 .75]);
scatter3(x,y,z,'MarkerEdgeAlpha',.8);
% [xData, yData, zData] = prepareSurfaceData(x,y,z);
% xyz=[xData yData zData]';
xyz=[x y z]';
% xyz=sortrows([x y z],[1 2])';
pieces = 15;
s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
pp1 = splinefit(s,xyz,pieces,0.75);
xyz1 = ppval(pp1,s);
hold on
plot3(xyz1(1,:),xyz1(2,:),xyz1(3,:),'-','linewidth',5);
% scatter3(xyz1(1,:),xyz1(2,:),xyz1(3,:));
grid on
% hold off
if nargout>0
    D=pdist2(xyz',xyz1');
    d=min(D,[],2);
end

