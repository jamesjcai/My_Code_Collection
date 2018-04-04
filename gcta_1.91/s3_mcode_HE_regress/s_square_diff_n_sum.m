function [Xsd,Xss,Xcp]=s_square_diff_n_sum(y)
n=length(y);
Xsd=zeros(n);
Xss=zeros(n);
for i=1:n-1
    for j=i+1:n
        %X(i,j)=(y(i)-y(j)).^2;
        Xsd(i,j)=(y(i)-y(j))*(y(i)-y(j));
        Xss(i,j)=(y(i)+y(j))*(y(i)+y(j));
    end
end
if nargout>2, 
    if size(y,1)==1
        Xcp=y'*y; 
    else
        Xcp=y*y'; 
    end
end
