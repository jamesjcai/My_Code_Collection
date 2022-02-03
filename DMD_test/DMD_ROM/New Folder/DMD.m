function [Phi,D] = DMD( X1 , X2 , n )

% X2 = A * X1

[U,S,V] = svd( X1,'econ' );

UU = U(:,1:n);
VV = V(:,1:n);
SS = S(1:n,1:n);
iSS = inv(SS);

AA = conj(UU') * X2 * VV * iSS;

[W,D] = eig( AA );

Phi = X2 * VV * iSS * W;