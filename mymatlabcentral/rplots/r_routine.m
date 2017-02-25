saveR('Data.R','x','y')
system([getRscript, ' ', mfilename('fullpath'), '.r']);
imshow('output.tif')