function [W,H]=nnmf_another(X,k)
if nargin<2, k=2; end

% https://www.mathworks.com/matlabcentral/answers/316708-non-negative-matrix-factorization
[m,n]=size(X); % n is # of samples, m is # of features
w0=rand(m,k);
h0=rand(k,n);
err=[];

%[r,c]=size(X); % c is # of samples, r is # of features
%H=rand(k,c);

maxIter=100;

for t=1:maxIter
    numer=w0'*X;
    H=h0.*(numer ./((w0'*w0)*h0 + eps(numer)));
    numer = X*h0';
    W=w0.*(numer ./(w0*(h0*h0') + eps(numer)));
    w0o=W; h0=H;
    err=[err norm(X-w0*h0)];
end
err
        