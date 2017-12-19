%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, a, dgp] = l0regdual(Q, z,  lam,  p, itn,  epn, eta, theta0)            
% Funnction [w, a, dgp] = lpregem(Q, z, lam, p, epn, eta, theta0)          
% takes the input and estimate theta0, a, and dgp respectively.
% Q: the input
% z: output
% theta0: initial theta0
% p: p in (0, 1]
% epn: small value
% itn: maximum interation number
% eta: sparsity control
%output: theta: premary parameter, a: dual parameters, dgp: gape between dual
%and primal solution.
%Zhenqiu Liu: Cedars-Sinai Medical Center
% 04/06/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, m] = size(Q);

if nargin < 8,
    theta0 = ones(m, 1);
end

if nargin < 7,
    eta = 1e-3;
end;

if nargin < 6,
    epn =1e-6;
end

if nargin < 5,
    itn = 1000;
end

if nargin < 4,
   p = 0;
end

theta = theta0;
I = diag(ones(n,1));
J = 1;
cont = 1;
%EM algorithm begins
while cont,
   u = theta;   % E-step
   % M-step:
   up = abs(u').^(2-p);
   Qu = Q.*repmat(up, n, 1);
   Ku = Q*Qu';
   a = (Ku +lam*I)\z;
   theta = Qu'*a;
   % theta(theta<0) =0;
   J = J +1;
   if J > itn,
       disp('Warning: No enough iterations');
       break;
   end
   if norm(theta-u) < epn,
       cont = 0;
   end;
   
end

dgp = norm(Ku*a -Q*theta);
theta(abs(theta) < eta) = 0;

   