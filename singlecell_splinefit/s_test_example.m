load E-MTAB-5989_dermal_fibroblasts.mat
%load GSM3204304_lung_airway_epithelial_cells.mat

[lgu0,dropr0,lgcv0,glist0]=i_scmatrix_stat(X0,genelist);
figure;
[d0,xyz0]=i_3dspline(lgu0,dropr0,lgcv0);

[lgu1,dropr1,lgcv1,glist1]=i_scmatrix_stat(X1,genelist);
[d1,xyz1]=i_3dspline(lgu1,dropr1,lgcv1);

%{
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
%}

[glist01,i,j]=intersect(glist0,glist1,'stable');
d0=d0(i); d1=d1(j); 
dd=d1-d0;
Tres2=table(glist01,d1,d0,dd);
Tres2=sortrows(Tres2,4,'descend');


[glistfinal,idx1,idx2]=intersect(glist01,genelist,'stable');
X0=X0(idx2,:);
X1=X1(idx2,:);


figure;
i_scatter_withinfo(d0,d1,glist01,X0,X1)
hline = refline(1);
hline.Color = 'r';


%% 
%{
figure;
i=d0<quantile(d0,0.1); 
j=d0>quantile(d0,0.90);
% scatter3(lgu0(i),dropr0(i),lgcv0(i),46,'g','markeredgealpha',0.5); 
scatter3(lgu0(j),dropr0(j),lgcv0(j),46,'b','markeredgealpha',0.5);
hold on
plot3(xyz1(1,:)',xyz1(2,:)',xyz1(3,:)','r-','linewidth',6)
grid on
% view(15,13)
view(70,10);
xlabel('log10(u)');
ylabel('Dropout Rate');
zlabel('log10(CV)');
%}