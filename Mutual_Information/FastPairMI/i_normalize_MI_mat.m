function [D]=i_normalize_MI_mat(D,methodid)
if nargin<2, methodid=2; end
n=size(D,1);
switch methodid
    case 1        
        d=zeros(n);
        for i=1:n-1
            for j=i+1:n
                d(i,j)=D(i,j)./max([D(i,i),D(j,j)]);
            end
        end
        D=d;
    case 2
        c=D./diag(D);
        c1=triu(c);
        c2=tril(c)';
        c(c1>c2)=c2(c1>c2);
        % c(1:n+1:end)=0;
        D=triu(c,1);
end