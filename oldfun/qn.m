function [x] = qn(x)

[x,y] = sort(x);
[~,z] = sort(y);
x = sort(mean(x,2));
x = x(z);
