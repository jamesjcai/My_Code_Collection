function [condI] = condI_estimation(Y,ds,co)
%function [condI] = condI_estimation(Y,ds,co)
%Estimates conditional mutual information I(y^1,...,y^M|y^{M+1}), where the m^th subspace is ds(m)-dimensional, using the specified method.
%
%INPUT:
%   Y: Y(:,t) is the t^th sample.
%  ds: subspace dimensions. ds(m) = dimension of the m^th subspace, m=1,...,M+1 (M+1=length(ds)), the last block is the conditioning variable.
%  co: conditional mutual information estimator object (structure).

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
 
%Here, we make use of the function handle initialized in 'condI_initialization':
    condI = co.function_handle(Y,ds,co); 