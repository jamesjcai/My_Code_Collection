clear all;

%Quick Test of DMD with control.

%Data collection

A = [1.5 0;0 0.1];
x0 = [4;7];
K = [-1];
m = 20;
DataX = x0;
DataU = [0];

B = [1;0];

for j = 1 : m 
   
  DataX(:,j+1) = A * (DataX(:,j)) + B.*(K*DataX(:,j));
  DataU(:,j) = K .* DataX(1,j);
end

%Data matrices
X   = DataX(:,1:end-1);
Xp  = DataX(:,2:end);
Ups = DataU;

%SVD
[U,Sig,V] = svd(X,'econ');

thresh = 1e-10;
r = length(find(diag(Sig)>thresh));

U    = U(:,1:r);
Sig  = Sig(1:r,1:r);
V    = V(:,1:r);

%DMD

A_DMD  = Xp*V*inv(Sig)*U'

%DMDc - B is known 

A_DMDc = (Xp - B*Ups)*V*inv(Sig)*U'
