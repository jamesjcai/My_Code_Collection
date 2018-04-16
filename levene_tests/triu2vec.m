function [x]=triu2vec(X,k)
% V = triu(A,k) returns the elements on and above the kth diagonal of A.
if nargin<2, k=1; end
mask=triu(true(size(X)),k);
x=X(mask);
% ref: https://stackoverflow.com/questions/29928468/transfer-the-lower-triangular-part-of-a-matrix-into-a-vector-in-matlab

