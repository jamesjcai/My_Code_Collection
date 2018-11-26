function r_heatmapcountour(X,G)
if nargin<2
    X=randn(3000,1)+10;
    Y=X.*1.5+randn(3000,1);
end
saveR('Data.R','X','Y')
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')


