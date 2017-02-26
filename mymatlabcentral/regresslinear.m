function [b1,b0,p,r2]=regresslinear(x,y)

if nargout<3
   [z]=polyfit(x,y,1); b1=z(1); b0=z(2);
else
   [z, ~, ~, ~, stats]=regress(y,[ones(size(x)) x]);
   b1=z(2);
   b0=z(1);
   p=stats(3);
   r2=stats(1);
end

