% ccacov.m - CCA applied to a covariance matrix
%
% Usage: [A,C,B] = ccacov(R,p)
%
% R = covariance matrix of size (p+q)x(p+q)
% p = size of Raa, the upper-left block of R 
%
% A,B = CCA coefficient matrices of sizes p x p and q x q
% C   = diagonal matrix of canonical correlations of size p x q
%
% notes: [Ua,Da,Va] = svd(Raa); Sa = sqrt(Da);  % eigenvalue decomposition of Raa
%        [Ub,Db,Vb] = svd(Rbb); Sb = sqrt(Db);
%        Cab = Sa \ Va' * Rab * Vb / Sb;
%        [Fa,C,Fb]  = svd(Cab);
%        A = Va / Sa * Fa = p x p matrix of coefficients
%        B = Vb / Sb * Fb = q x q matrix of coefficients
%        
%        y = [ya; yb] = (p+q)x1 dimensional vector
%        R = E[y^* y^T] = [Raa Rab
%                          Rba Rbb]
%
%        canonical basis: 
%        wa = A^T ya, wb = B^T yb,  w = [wa; wb], 
%        Rww = E[w^* w^T] = [Ip  C 
%                            C'  Iq]
%
%        A' * Raa * A = Ip = p x p identity matrix
%        A' * Rab * B = C
%        B' * Rbb * B = Iq
%
%        see also CCA - applied to data matrices Ya,Yb
%
%        with Y=[Ya,Yb] and R=Y'*Y, ccacov(R,p) gives same result as cca(Ya,Yb)

% S. J. Orfanidis - 1999
% ECE Department
% Rutgers University
% Piscataway, NJ 08854
% email: orfanidi@ece.rutgers.edu

function [A,C,B] = ccacov(R,p)

if nargin==0, help ccacov; return; end

M = size(R,1);  

Raa = R(1:p, 1:p);          % upper-left  p x p sub-block
Rab = R(1:p, p+1:M);        % upper-right p x q sub-block
Rbb = R(p+1:M, p+1:M);      % lower-right q x q sub-block

[Ua,Da,Va] = svd(Raa);      % eigenvalue decomposition of Raa, Ua=Va
[Ub,Db,Vb] = svd(Rbb);

Sa = sqrt(abs(Da));         % theoretically, Da,Db are positive
Sb = sqrt(abs(Db));

Cab = Sa \ Va' * Rab * Vb / Sb;

[Fa,C,Fb]  = svd(Cab);

A = Va / Sa * Fa;
B = Vb / Sb * Fb;


