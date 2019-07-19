function [T]=sc_3dfeadev2x(X,Y,genelist)

i=sum(X,2)==0 & sum(Y,2)==0;
if sum(i)>0
    X(i,:)=[]; Y(i,:)=[]; genelist(i)=[];
    warning('%d lowly-expressed genes removed.',sum(i));
end
% i=find(sum(X,2)==0);
% j=randi([1 size(X,2)],length(i),1);
% X(sub2ind(size(X),i,j))=1;
% 
% i=find(sum(Y,2)==0);
% j=randi([1 size(Y,2)],length(i),1);
% Y(sub2ind(size(Y),i,j))=1;

[T0,pp0]=sc_3dfeadev(X,genelist);
glist0=T0.genes;
d0=T0.d;

[x,y,z,glist1]=sc_3dfeastat(Y,genelist);
xyz=[x y z]';
s = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)]);
xyz1 = ppval(pp0,s);
D=pdist2(xyz',xyz1');
d1=min(D,[],2);


%[T1]=sc_3dfeadev(Y,genelist);
%d0=T0.d; glist0=T0.genes;
%d1=T1.d; glist1=T1.genes;
%
[glist01,i,j]=intersect(glist0,glist1,'stable');
d0=d0(i); d1=d1(j); 
dd=d1-d0;
genes=glist01;
T=table(genes,d0,d1,dd);
T=sortrows(T,'dd','descend');



