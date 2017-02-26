function r_densityhist2(x1,y1,x2,y2)
if nargin<1
    % Generate some random data
    x1=5+randn(1,1000)*2;
    y1=randn(1,1000);
    x2=3+randn(1,1000)*2;
    y2=2+randn(1,1000);    
end
saveR('Data.R','x1','y1','x2','y2')
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')
