function B = l0regnet(X1,X2, epn, itn, choice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X1 is an NxM matrix, N is the number of samples
%M is the number of Glycans
%X2, is an NxM1 matrix, and M1 is the number of SNPs
% Y: Nx1  class vector (0/1)
% lambda is penalty for lasso. It can be a vector
%B is the matrix of PxP with nonzero diagnal if it is differentiated
%E: is the minimal MSE for all genes (px1)
%oplamb: Px1 vector for the optimal values
% Written by Zhenqiu Liu
% Cedars Sinai Medical Center
% 07/02/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
if size(X1,1) < 2
    error(message('L0:TooFewObservations'));
end

if size(X1, 1) ~= size(X2, 1),
   error(message(' The sample  size of Glycan and SNP must be equal'));
end


% If Y is a row vector, convert to a column vector

if nargin < 5,
    choice = 1;
end

if nargin < 4,
    itn = 1000;
end

if nargin < 3,
    epn = 1e-6;
end

m1 = size(X1,2);
m2 = size(X2, 2);
m = m1+ m2;
B =zeros(m, m);

  for j = 1:m1,
      %j 
      Zin = X2;
      Zo =X1(:,j);
  
      % Standardize input and output
     [n, m] = size(Zin);
     %Zo = Zo - mean(Zo);
     %Zin = (Zin - repmat(mean(Zin), n1, 1))./repmat(std(Zin), n1, 1);
      if choice ==1;  %AIC
          lam = 2;
      elseif choice ==2,  %BIC
          lam = log(n);
      else
          lam = log(m);  % RIC Risk Inflation Criteria
      end
      if n <= m 
         b = l0regdual(Zin, Zo,  lam,  0, itn, epn);
      else
          b = l0regprimal0(Zin, Zo,  lam,  0, itn, epn);
      end
      B(:,j) = [zeros(m1, 1); b];
      j
  end   
  B = (B + B');
  k
end

