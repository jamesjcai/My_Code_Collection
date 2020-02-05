function [b]=i_head(filename,n)
if nargin<2, n=50; end
fid=fopen(filename,'r');
b=textscan(fid,'%s',n,'delimiter','\n');
fclose(fid);
b=string(b{1});