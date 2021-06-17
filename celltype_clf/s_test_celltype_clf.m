T=readtable('G1_singlecells_counts.txt');
x1=table2array(T(:,5:end));
T=readtable('G2M_singlecells_counts.txt');
x2=table2array(T(:,5:end));
T=readtable('S_singlecells_counts.txt');	
x3=table2array(T(:,5:end));
g=string(T.AssociatedGeneName);
idx=g~="NA";
X=[x1 x2 x3];
X=X(idx,:);
g=g(idx);
id=ones(96,1);
idv=[id;2*id;3*id];

%%
load('data_celltype_clf.mat')
[X,g]=sc_selectg(X,g,50);
[T,X,g]=sc_hvg(X,g);

X=sc_transform(X);
Xt=X(1:2000,:)';








