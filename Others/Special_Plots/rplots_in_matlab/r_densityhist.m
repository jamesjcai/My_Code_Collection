function r_densityhist(x)
if nargin<1
    % Generate some random data
    x=randn(1,1000);
end
saveR('Data.R','x')
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')
