function Phi=Gaussian_kernel(Phi_tmp,sigma)

Phi = exp(Phi_tmp/sigma^2);
