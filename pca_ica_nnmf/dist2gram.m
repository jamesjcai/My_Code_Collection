function [G]=dist2gram(D)
n=size(D,1);
v1=ones(n,1);
D=D.*D;
J=eye(n)-(1/n)*(v1*v1');
G=-0.5*(J*D*J);

