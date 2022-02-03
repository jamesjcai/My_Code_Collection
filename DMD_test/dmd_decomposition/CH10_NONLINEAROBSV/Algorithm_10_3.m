X1 = X(:,1:end-1); %% data snapshots
X2 = X(:,2:end); %% shifted data

[U2,Sigma2,V2] = svd(X1, 'econ');
r=10;  %% rank 10 truncation
U=U2(:,1:r); Sigma=Sigma2(1:r,1:r); V=V2(:,1:r);
Atilde = U'*X2*V/Sigma;
[W,D] = eig(Atilde);
Phi=X2*V/Sigma*W;

lambda=diag(D);
omega=log(lambda)/dt;

y0 = Phi\X(:,1);  %% project on initial conditions

q_modes = zeros(r,length(t));  %% DMD modes
for iter = 1:length(tf)
    q_modes(:,iter) =(y0.*exp(omega*(tf(iter))));
end

q_dmd = Phi*q_modes;   %% DMD approximation