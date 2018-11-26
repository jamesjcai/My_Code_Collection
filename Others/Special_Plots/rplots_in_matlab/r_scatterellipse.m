function r_scatterellipse(x,y)
if nargin<2
    x=randn(1000,1);
    y=randn(1000,1)+0.5*rand(1000,1)+x;
end

%%
saveR('Data.R','x','y');
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')

