function [v]=frechetvar(X,u0)
v=sum(vecnorm(X-u0,2,2));