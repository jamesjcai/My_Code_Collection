function r_linepointerrorbar(xx,yy)
if nargin<2
    xx=[1:12 1:12]';
    yy=[8, 12, 13, 18,  22, 16, 24, 29,  34, 15, 8, 6,...
                         9, 10, 12, 18, 26, 28, 28, 30, 20, 10, 9, 9]';
end

%%
saveR('Data.R','xx','yy')
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')

