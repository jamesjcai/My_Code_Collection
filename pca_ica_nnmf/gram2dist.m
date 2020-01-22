function [D]=gram2dist(G)
n=size(G,1);
v1=ones(n,1);
H=v1*v1'*diag(diag(G));
D=H+H'-2*G;
D=sqrt(D);

