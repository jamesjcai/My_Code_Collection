function D = interdist(A,methodid)
if nargin<2, methodid=1; end
switch methodid
    case 1
        n=size(A,2);
        for k=1:n
            D(:,:,k)=A(:,k)-A(:,k)';
        end
        D=sqrt(sum(D.^2,3));
    case 2
        % D = sqrt(bsxfun(@plus,sum(A.^2)',sum(A.^2))-2*(A'*A));
        % https://rundong.wordpress.com/2013/09/09/efficient-matlab-i-pairwise-distances/
        D=sqrt(sum(A.^2,2)+sum(A.^2,2)'-2*(A*A'));
end
end
