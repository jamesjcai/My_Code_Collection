load zxy Zc X Y
Z=cellfun(@(X)(sum(X<0.001)),Zc);
figure; mesh(X,Y,Z)
figure; contour(X,Y,Z,'showtext','on')
% https://stackoverflow.com/questions/44816434/matlab-how-to-make-smooth-contour-plot
