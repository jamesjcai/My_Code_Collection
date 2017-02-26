function c=equsumpooling(a,k)
%a=[3 7 8 2 5 7 6 4 1 9];
c=zeros(size(a));
m=sum(a)/k;
x=cumsum(a);
for (j=1:k)
 c=c+j*(x>(j-1)*m & x<=j*m);
end
