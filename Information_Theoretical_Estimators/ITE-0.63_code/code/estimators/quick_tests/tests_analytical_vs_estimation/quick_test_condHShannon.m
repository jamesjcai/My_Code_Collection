%function [] = quick_test_condHShannon()
%Quick test for conditional Shannon entropy estimators: analytical expression vs estimated value as a function of the sample number. In the test, normal variables are considered. See also 'analytical_value_condHShannon.m'.

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

%clear start:
    clear all; close all;

%parameters:
    distr = 'normal'; %fixed
    d1 = 1; %dimension of the distribution-1
    d2 = 2; %dimension of the distribution-2
    num_of_samples_v = [1000:1000:30*1000]; %sample numbers used for estimation
    %estimator:
        %meta:
            cost_name = 'Shannon_HShannon';  %dm>=1
    
%initialization:    
    num_of_samples_max = num_of_samples_v(end);
    L = length(num_of_samples_v);
    co = condH_initialization(cost_name,1);
    condH_hat_v = zeros(L,1); %vector of estimated conditional entropies
    d = d1 + d2;
    par.d1 = d1;
    par.d2 = d2;
    
%distr, d -> samples (Y=[Y1;Y2]), analytical value [condH = H(y^1|y^2)]:
    switch distr 
        case 'normal'
            %mean:
                m = rand(d,1);
            %random linear transformation applied to N(0,I):
                A = rand(d); 
                %A = eye(d); %do not transform the data
            %covariance matrix:
                C = A * A.';
            %generate samples:
                Y = A * randn(d,num_of_samples_max) + repmat(m,1,num_of_samples_max); %AxN(0,I)+m
            par.cov = C;
        otherwise
            error('Distribution=?');
    end  
    %analytical value:
        condH = analytical_value_condHShannon(distr,par);
            
%estimation:
    Tk = 0;%index of the sample number examined
    for num_of_samples = num_of_samples_v
        Tk = Tk + 1;
        condH_hat_v(Tk) = condH_estimation(Y(1:d1,1:num_of_samples),Y(d1+1:d,1:num_of_samples),co);
        disp(strcat('Tk=',num2str(Tk),'/',num2str(L)));
    end
    
%plot:
    plot(num_of_samples_v,condH_hat_v,'r',num_of_samples_v,condH*ones(L,1),'g');
    legend({'estimation','analytical value'});
    xlabel('Number of samples');
    ylabel('Conditional Shannon entropy');
