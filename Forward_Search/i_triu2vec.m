function [x]=i_triu2vec(X,k)
if nargin<2, k=0; end
x=X(triu(true(size(X)),k));