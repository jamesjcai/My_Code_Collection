% https://www.nature.com/articles/nmeth.2645#s1
% page 24.
load examgrades.mat
% x=X(1:10,1:20);
% gm2 = @(x) exp(sum(log(x))/numel(x));
% sf=nanmedian(x./nangeomean(x,2));
% sf1=median(x./geomean(x,2));
% sf2=median(x./geomean(x,2));
% x_n=x./sf;
x=X;
        xn0=x;
        xn0(xn0==0)=nan;
        sf=nanmedian(xn0./nangeomean(xn0,2));
        x=x./sf;