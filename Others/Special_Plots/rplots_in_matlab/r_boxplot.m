function r_boxplot(X,G)
if nargin<2
    X=[randn(100,1); 0.2+randn(100,1)];
    G=[ones(100,1); 2*ones(100,1)];
end

saveR('Data.R','X','G')
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')