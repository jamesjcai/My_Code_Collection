xtilde = zeros(n,n);
XDAT = [];
XDATtilde = [];
XDATtildeNoise = [];
for t = 0:dt:T
    xtilde = 0*xtilde;
    for k=1:K
        xtilde(I(k),J(k)) = exp(damping(k)*t)*cos(2*pi*F(k)*t)*IC(k) + exp(damping(k)*t)*sqrt(-1)*sin(2*pi*F(k)*t)*IC(k);
    end
    XDATtilde = [XDATtilde reshape(xtilde,n^2,1)];
    xRMS = sqrt((1/(n*n))*sum(xtilde(:).^2));
    xtilde = xtilde + noisemag*xRMS*randn(size(xtilde)) + noisemag*xRMS*sqrt(-1)*randn(size(xtilde));
    XDATtildeNoise = [XDATtildeNoise reshape(xtilde,n^2,1)];
    x = real(ifft2(xtilde));
    XDAT = [XDAT reshape(x,n^2,1)];
end
if(saveFLAG)
    save([filename,'_DATAX.mat']);
end