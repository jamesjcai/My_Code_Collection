function Phi=Epanechnikov_kernel(Phi_tmp,sigma)

invsigma = 1/sigma^2;
Phi = 1 + max(-1,Phi_tmp*invsigma);
