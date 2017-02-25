function [s]=fillsfs(x,n,value,binvalue)

s=zeros(1,n);

if nargin==2
    for k=1:n
        s(k)=sum(x==k);    
    end
elseif nargin==3
    for k=1:n
        s(k)=mean(value(x==k));
    end
elseif nargin==4
    s=zeros(1,length(binvalue));    
    for k=1:length(binvalue)
        s(k)=sum(x==binvalue(k));
    end    
end