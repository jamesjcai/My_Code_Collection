load zxy Zc X Y
Z=cellfun(@(X)(sum(X<0.001)),Zc);
figure;
mesh(X,Y,Z)
