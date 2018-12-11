load DermalFibroblasts_data.mat

[lgu0,dropr0,lgcv0,glist0]=i_scmatrix_stat(data0,glist);
d0=i_3dspline(lgu0,dropr0,lgcv0);

[lgu1,dropr1,lgcv1,glist1]=i_scmatrix_stat(data1,glist);
d1=i_3dspline(lgu1,dropr1,lgcv1);

[glist01,i,j]=intersect(glist0,glist1,'stable');
d0=d0(i); d1=d1(j); 
dd=d1-d0;
Tres=table(glist01,d1,d0,dd);
Tres=sortrows(Tres,4,'descend');

figure;
i_scatter_withinfo(d0,d1,glist01)