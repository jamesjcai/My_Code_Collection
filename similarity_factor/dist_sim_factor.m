function [s,d]=dist_sim_factor(Xs,Xh)

u=mean(Xh)-mean(Xs);
[sig]=robustcov(Xs);
d=u*sig*u';

fun = @(z)exp(-z.^2/2);
m=sqrt(2/pi);
s = m*integral(fun,d,Inf);