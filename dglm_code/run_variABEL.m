function [p,effsize]=run_variABEL(y1,x)

    [b]=regress(y1,[ones(size(x)) x]);
    y2=y1-(b(1)+b(2)*x);
    [effs, ~, ~, ~, Stats]=regress(y2.^2,[ones(size(x)) x]);
    p=Stats(3);
    effsize=effs(2);