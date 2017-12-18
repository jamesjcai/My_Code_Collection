function [condH] = condHShannon_HShannon_estimation(Y1,Y2,co)
%function [condH] = condHShannon_HShannon_estimation(Y1,Y2,co)
%Estimates the conditional Shannon entropy of Y1 given Y2 [H(y^1|y^2)] via the relation H(y^1|y^2) = H([y^1;y^2]) - H(y^2) relation, where H is the Shannon differential entropy.
%
%Note:
%   1)We use the naming convention 'condH<name>_estimation' to ease embedding new conditional entropy estimation methods.
%   2)This is a meta method: the Shannon entropy estimator can be arbitrary.
%
%INPUT:
%   Y1: Y1(:,t1) is the t1^th sample.
%   Y2: Y2(:,t2) is the t2^th sample.
%  co: entropy estimator object.

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

num_of_samplesY1 = size(Y1,2);
num_of_samplesY2 = size(Y2,2);
    
%Shannon entropy of y^2:    
    H2 =  H_estimation(Y2,co.member_co); 
    
%Shannon entropy of [y^1;y^2]:
    num_of_samples = min(num_of_samplesY1,num_of_samplesY2);
    H12 = H_estimation([Y1(:,1:num_of_samples);Y2(:,1:num_of_samples)],co.member_co); 
    
condH = H12 - H2;

