clear all;
close all;

X   = StateData(:,1:end-1);
Xp  = StateData(:,2:end);
Ups = InputData(:,1:end-1);

Omega = [X;Ups];

[U,Sig,V] = svd(Omega,'econ');

thresh = 1e-10;
rtil = length(find(diag(Sig)>thresh));

Util    = U(:,1:rtil); 
Sigtil  = Sig(1:rtil,1:rtil);
Vtil    = V(:,1:rtil); 

[U,Sig,V] = svd(Xp,'econ');

thresh = 1e-10;
r = length(find(diag(Sig)>thresh));

Uhat    = U(:,1:r); 
Sighat  = Sig(1:r,1:r);
Vbar    = V(:,1:r); 

n = size(X,1); 
q = size(Ups,1);
U_1 = Util(1:n,:);
U_2 = Util(n+q:n+q,:);

approxA = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_1'*Uhat;
approxB = Uhat'*(Xp)*Vtil*inv(Sigtil)*U_2';

[W,D] = eig(approxA);

Phi = Xp * Vtil * inv(Sigtil) * U_1'*Uhat * W;