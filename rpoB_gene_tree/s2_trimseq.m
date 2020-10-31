x=fastaread('aaa.fas');
for k=2:length(x)
    k
[a,b,c]=swalign(x(1).Sequence,x(2).Sequence);
disp(c(1))
x(k).Sequence=upper(x(k).Sequence(c(2)-10:c(2)+500));
end
fastawrite('bbbbb2.fas',x);
