function [g]=split2bins(x,n)

l=length(x);
g=zeros(l,1);
offset=rand(l,1).*range(x)./1000;
x=x+offset;
borders=prctile(x,[0:(100/n):100]);


for k=1:n
    g(x>=borders(k)&x<borders(k+1))=k;    
end






