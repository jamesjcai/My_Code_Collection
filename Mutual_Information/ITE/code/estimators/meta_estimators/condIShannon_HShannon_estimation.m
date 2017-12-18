function [condI] = condIShannon_HShannon_estimation(Y,ds,co)
%function [condI] = condIShannon_HShannon_estimation(Y,ds,co)
%Estimates conditional Shannon mutual information (condI) using Shannon differential entropy and 
%the I(y^1,...,y^M|y^{M+1}) = -H([y^1;...;y^M;y^{M+1}]) + \sum_{m=1}^M H([y^m;y^{M+1}]) - (M-1) H(y^{M+1}) relation.
%
%Note:
%   1)We use the naming convention 'condI<name>_estimation' to ease embedding new conditional mutual information estimation methods.
%   2)This is a meta method: the Shannon differential entropy estimator can be arbitrary. 
%
%INPUT:
%   Y: Y(:,t) is the t^th sample.
%  ds: subspace dimensions. ds(m) = dimension of the m^th subspace, m=1,...,M+1 (M+1=length(ds)), the last block is the conditioning variable.
%  co: conditional mutual information estimator object.

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

%co.mult:OK. The information theoretical quantity of interest can be (and is!) estimated exactly [co.mult=1]; the computational complexity of the estimation is essentially the same as that of the 'up to multiplicative constant' case [co.mult=0]. In other words, the estimation is carried out 'exactly' (instead of up to 'proportionality').

%verification:
    if sum(ds) ~= size(Y,1);
        error('The subspace dimensions are not compatible with Y.');
    end
    if length(ds) <= 2;
        error('At least two non-conditioning subspaces are needed.');
    end
     
%cum_ds, M:
    cum_ds = cumsum([1;ds(1:end-1)]);%1,d_1+1,d_1+d_2+1,...,d_1+...+d_{M}+1 = starting indices of the subspaces (M+1=number of subspaces).
    M = length(ds)-1;

%H_joint:    
    H_joint = H_estimation(Y,co.member_co);
    
%H_cross:    
    H_cross = 0;
    indMp1 = [cum_ds(M+1):cum_ds(M+1)+ds(M+1)-1];
    for m = 1 : M
        indm = [cum_ds(m):cum_ds(m)+ds(m)-1];
        H_cross = H_cross + H_estimation([Y(indm,:);Y(indMp1,:)],co.member_co);
    end
    
%H_condition:    
    H_condition = H_estimation(Y(indMp1,:),co.member_co);
    
condI = -H_joint + H_cross - (M-1) * H_condition;
