function [x]=s_extract_triu(X)
n=size(X,1);
x=[];
for k=1:n-1
    x=[x;diag(X,k)];
end