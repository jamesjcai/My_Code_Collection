function [y]=ls2cellstr(dirname,isdir)
if nargin<2
    isdir=true;    
end
x=ls(dirname);
if isempty(x)
    y='';
end
n=size(x,1);

if isdir
    for k=3:n
        y{k-2}=strtrim(x(k,:));
    end
else
    for k=1:n
        y{k}=strtrim(x(k,:));
    end    
end