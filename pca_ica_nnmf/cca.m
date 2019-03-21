% cca.m - Canonical Correlation Analysis
%
% Usage: [A,C,B] = cca(Ya,Yb)
%
% Ya,Yb = data matrices of sizes N x p and N x q
%
% A,B = CCA coefficient matrices of sizes p x p and q x q
% C   = diagonal matrix of canonical correlations of size p x q
%
% notes: [Ua,Sa,Va] = svd(Ya, 0)
%        [Ub,Sb,Vb] = svd(Yb, 0)
%        [Fa,C,Fb]  = svd(Ua'*Ub)
%        A = Va / Sa * Fa = p x p matrix of coefficients
%        B = Vb / Sb * Fb = q x q matrix of coefficients
%
%        Wa = Ya*A = N x p with orthonormal columns, Wa'*Wa = Ip
%        Wb = Yb*B = N x q with orthonormal columns, Wb'*Wb = Iq
%        Wa'*Wb = C = diagonal matrix of canonical correlation coefficients
%
%        [Wa,Wb] has correlation matrix: [Wa,Wb]'*[Wa,Wb] = [Ip, C; C', Iq] 
%
%        Ya,Yb must have zero-mean columns (use ZMEAN)
%     
%        the canonical angles between the subspaces Ya,Yb are theta=acos(diag(C))
%        the maximum canonical angle corresponding to minimum C is the same as from 
%        the built-in MATLAB function SUBSPACE(Ya,Yb)
%
%        see also CCACOV - applied to a given covariance matrix
%
%        with Y=[Ya,Yb] and R=Y'*Y, ccacov(R,p) gives same result as cca(Ya,Yb)

% S. J. Orfanidis - 1999
% ECE Department
% Rutgers University
% Piscataway, NJ 08854
% email: orfanidi@ece.rutgers.edu

function [A,C,B] = cca(Ya,Yb)

if nargin==0, help cca; return; end

[Ua, Sa, Va] = svd(Ya, 0);
[Ub, Sb, Vb] = svd(Yb, 0);
[Fa, C,  Fb]  = svd(Ua'*Ub);
A = Va / Sa * Fa;
B = Vb / Sb * Fb;


