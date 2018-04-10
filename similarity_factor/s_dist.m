function [ s3,s4 ] = s_dist(Xs,Xh)
xh=mean(Xh);
xs=mean(Xs);

%%Calculation of seudo-inverse of Snapshot
A=Xs;
 [U,S,V] = svd(A); % A = U*S*V'
 T=S;
 T(find(S~=0)) = 1./S(find(S~=0));
xfc=T(1:min(size(Xs)),1:min(size(Xs)));


s3=sqrt((xh-xs)*xfc*(xh-xs)');% Equation 3

syms x;
s=sqrt(2/pi)*int(exp(-x.^2/2),x,s3,inf);
s4=vpa(s,4);%%%Equation 4

end

