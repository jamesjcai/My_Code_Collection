function [T]=sc_3dfeadev2(X,Y,genelist,showit)
if nargin<4, showit=false; end

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

[T0]=sc_3dfeadev(X,genelist);
[T1]=sc_3dfeadev(Y,genelist);

d0=T0.d; glist0=T0.genes;
d1=T1.d; glist1=T1.genes;

[glist01,i,j]=intersect(glist0,glist1,'stable');
d0=d0(i); d1=d1(j); 
dd=d1-d0;
genes=glist01;
T=table(genes,d0,d1,dd);
T=sortrows(T,'dd','descend');


if showit
    [~,~,idx]=intersect(T.genes,genelist,'stable');
    X0trim=X(idx,:);
    X1trim=Y(idx,:);
    i_scatter_withinfo(T.d0,T.d1,T.genes,X0trim,X1trim)
    hline = refline(1);
    hline.Color = 'r';
end
end

function i_scatter_withinfo(d0,d1,infotxt,X,Y)
    if nargin<5, Y=[]; end
    if nargin<4, X=[]; end
    scatter(log(d0+0.01),log(d1+0.01),[],log(d1)-log(d0));
    dt = datacursormode;
    if isempty(X)
        dt.UpdateFcn = {@i_myupdatefcn1,infotxt};
    else
        dt.UpdateFcn = {@i_myupdatefcn3,infotxt,X,Y};
    end
end


