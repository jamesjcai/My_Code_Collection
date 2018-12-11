load DermalFibroblasts_data.mat

[lgu0,dropr0,lgcv0,glist0]=i_scmatrix_stat(data0,glist);
figure;
d0=i_3dspline(lgu0,dropr0,lgcv0);

[lgu1,dropr1,lgcv1,glist1]=i_scmatrix_stat(data1,glist);
d1=i_3dspline(lgu1,dropr1,lgcv1);

[a, b]=view;
loops=length(10:5:200);
F(loops) = struct('cdata',[],'colormap',[]);
c=1;
v = VideoWriter('snakecylinder.avi');
open(v);
for k=10:5:200
    view(a+k,b);
    %drawnow
    F(c) = getframe(gcf);
    frame=getframe(gcf);
    writeVideo(v,frame);
    c=c+1;
    pause(0.5);
end
% figure
% movie(F,loops)
close(v)

[glist01,i,j]=intersect(glist0,glist1,'stable');
d0=d0(i); d1=d1(j); 
dd=d1-d0;
Tres=table(glist01,d1,d0,dd);
Tres=sortrows(Tres,4,'descend');

figure;
i_scatter_withinfo(d0,d1,glist01)
hline = refline(1);
hline.Color = 'r';
