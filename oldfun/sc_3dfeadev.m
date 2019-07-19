function [T,pp1]=sc_3dfeadev(X,genelist,sortit)
if nargin<3, sortit=true; end
[x,y,z,genes]=sc_stat(X,genelist);
lgu=x;
dropr=y;
lgcv=z;

xyz=[x y z]';
pieces = 15;
s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
pp1 = splinefit(s,xyz,pieces,0.75);
xyz1 = ppval(pp1,s);

D=pdist2(xyz',xyz1');
d=min(D,[],2);

T=table(genes,lgu,dropr,lgcv,d);
% 'variablenames',{'Genes','Log10_Mean','Dropout_Rate','Log10_CV','Deviation_3D'});
if sortit
    T=sortrows(T,'d','descend');
end