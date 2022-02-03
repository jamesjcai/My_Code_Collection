Y1=[X1; (X1.*abs(X1).^2)];  %% Data observables g1
Y2=[X2; (X2.*abs(X2).^2)];  %% Shifted Data

[U2,Sigma2,V2] = svd(Y1, 'econ');
r=10;  %% rank 10 truncation
U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
Atilde = U'*Y2*V/Sigma;
[W,D] = eig(Atilde);
Phi2=Y2*V/Sigma*W;