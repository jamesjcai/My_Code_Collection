function [condH] = condH_estimation(Y1,Y2,co)
%function [condH] = condH_estimation(Y1,Y2,co)
%Conditional entropy estimation (condH) of Y using the specified conditional entropy estimator.
%
%INPUT:
%   Y1: Y(:,t1) is the t1^th sample.
%   Y2: Y(:,t2) is the t2^th sample.
%  co: conditional entropy estimator object.

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

%Here, we make use of the function handle initialized in 'condH_initialization':
    condH = co.function_handle(Y1,Y2,co);