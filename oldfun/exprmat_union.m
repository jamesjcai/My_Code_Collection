function [X0o,X1o,g]=exprmat_union(X0,X1,g0,g1)
[g,i0,i1]=intersect(g0,g1,'stable');
m0=size(X0,2);
m1=size(X1,2);
n=size(g,1);
X0o=zeros(n,m0);
X1o=zeros(n,m1);
X0o(i0,:)=X0;
X1o(i1,:)=X1;






    