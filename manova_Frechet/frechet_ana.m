function [p] = frechet_ana(X1,X2)

n1=size(X1,1);
[~,u1]=kmeans(X1,1);
d1=sqrt(sum((X1-u1).^2,2));
v1=sum(d1.^2)/n1;
q1=sum(d1.^4)/n1-(v1).^2;

n2=size(X2,1);
[~,u2]=kmeans(X2,1);
d2=sqrt(sum((X2-u2).^2,2));
v2=sum(d2.^2)/n2;
q2=sum(d2.^4)/n2-(v2).^2;

% Pooled sample

X=[X1;X2];
n=size(X,1);
[~,u]=kmeans(X,1);
v=sum(sum((X-u).^2,2))/n;

%% Fn + Un
l1=n1/n;
l2=n2/n;

Fn=v-(l1*v1+l2*v2); %%%eq18
Un=l1*l2*(v1-v2).^(2)/(q1.*q2); %%%eq19

tn=n*Un/(l1/q1+l2/q2)+n*Fn.^2/(l1.^2*q1+l2.^2*q2); %%%eq22
p = chi2cdf(tn,1,'upper'); %%%eq23

