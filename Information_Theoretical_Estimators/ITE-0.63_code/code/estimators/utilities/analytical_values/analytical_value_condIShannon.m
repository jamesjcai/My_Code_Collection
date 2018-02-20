function [condI] = analytical_value_condIShannon(distr,par)
%function [condI] = analytical_value_condIShannon(distr,par)
%Analytical value of the conditional Shannon mutual information (condI) for the given distribution. See also 'quick_test_condIShannon.m'.
%
%INPUT:
%   distr  : name of the distribution.
%   par    : parameters of the distribution (structure).
%
%If distr = 'normal': par.ds = vector of component dimensions; par.cov = (joint) covariance matrix.

%Copyright (C) 2012- Zoltan Szabo ("http://www.gatsby.ucl.ac.uk/~szabo/", "zoltan (dot) szabo (at) gatsby (dot) ucl (dot) ac (dot) uk")
%
%This file is part of the ITE (Information Theoretical Estimators) toolbox.
%
%ITE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
%This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with ITE. If not, see <http://www.gnu.org/licenses/>.

%ds, cum_ds, M:
    ds = par.ds;
    cum_ds = cumsum([1;ds(1:end-1)]);%1,d_1+1,d_1+d_2+1,...,d_1+...+d_{M}+1 = starting indices of the subspaces (M+1=number of subspaces).
    M = length(ds)-1;
    
switch distr
    case 'normal'
        C = par.cov;        
        %H_joint:    
            H_joint = analytical_value_HShannon(distr,par);
        %H_cross:    
            H_cross = 0;
            indMp1 = [cum_ds(M+1):cum_ds(M+1)+ds(M+1)-1];
            for m = 1 : M
                indm = [cum_ds(m):cum_ds(m)+ds(m)-1];
                par.cov = C([indm,indMp1],[indm,indMp1]);
                H_cross = H_cross + analytical_value_HShannon(distr,par);
            end
        %H_condition:    
            par.cov = C(indMp1,indMp1);
            H_condition = analytical_value_HShannon(distr,par);
            
        condI = -H_joint + H_cross - (M-1) * H_condition;        
    otherwise
        error('Distribution=?');            
end