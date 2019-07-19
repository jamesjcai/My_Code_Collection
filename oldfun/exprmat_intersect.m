function [X0,X1,g]=exprmat_intersect(X0,X1,g0,g1)
[~,i0,i1]=intersect(g0(:,1),g1(:,1),'stable');
g=g0(i0,:);
X0=X0(i0,:);
X1=X1(i1,:);
    