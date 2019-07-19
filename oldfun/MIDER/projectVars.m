%--------------------------------------------------------------------------
% projectVars: project variables to a lower dimensional space
% 
% [Y] = projectVars(dissimilarities,ntotal) calculates a lower-dimensional 
% projection of the distances among variables, so that they are preserved 
% as much as possible, using a technique known as classical 
% Multidimensional Scaling (MDS) or Principal Coordinates Analysis (PCA). 
% The 2 first columns of the Y array are afterwards used as the horizontal 
% and vertical coordinates of the distance map (in the script plotResults.m)
%
% dissimilarities = ntotal*ntotal array of pairwise distances between variables
% ntotal          = number of variables
% Y               = ntotal*(ntotal-1) array of coordinates of the projected variables
%
% Dependencies: projectVars.m is called by plotResults.m
%
% Reference: 
% Seber, G. A. (2009). Multivariate observations. John Wiley & Sons.
%
% Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
% Created: 04/11/2014, last modified: 24/11/2014
%--------------------------------------------------------------------------
% Copyright (C) 2014  Alejandro Fernández Villaverde
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

function [Y] = projectVars(dissimilarities,ntotal)

% Reshape dissimilarities as array:
D = squareform(dissimilarities); 

% Calculate eigenvalues: 
W = eye(ntotal) - repmat(1/ntotal,ntotal,ntotal);
M = W * (-.5 * D .* D) * W;
[V E] = eig((M+M')./2);  

% Sort eigenvalues in descending order:
[e i] = sort(diag(E)); 
fe = flipud(e); 
fi = flipud(i); 

% Use only positive eigenvalues
pe = (fe > 0); 
Y = V(:, fi(pe)) * diag(sqrt(fe(pe)));  
    
end
