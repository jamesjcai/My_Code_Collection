function [T,W,H]=sc_nnmf(X,genelist,k)
if nargin<3, k=10; end
if nargin<2, genelist=[]; end
A=log(X+1);
[W,H]=nnmf(A,k);
if ~isempty(genelist)
    w=[];
    for wk=1:k
        [~,i]=sort(W(:,wk),'descend');
        w=[w genelist(i)];
    end
    T=array2table(w);
end
