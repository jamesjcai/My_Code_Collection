function [x]=i_extract_triu(X,methodid)

if nargin<2, methodid=3; end

switch methodid
    case 1
        warning('output vector element order.')
        n=size(X,1);
        x=[];
        for k=1:n-1
            x=[x;diag(X,k)];
        end
    case 2
        mask = bsxfun(@lt,[1:size(X,1)]',1:size(X,2));
        x=X(mask);
    case 3
        mask=triu(true(size(X)),1);
        x=X(mask);
end

% ref: https://stackoverflow.com/questions/29928468/transfer-the-lower-triangular-part-of-a-matrix-into-a-vector-in-matlab

