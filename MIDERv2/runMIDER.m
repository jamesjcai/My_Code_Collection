%--------------------------------------------------------------------------
% runMIDER: Call the MIDER algorithm (mider.m), save results, and plot them
% 
% Inputs: runMIDER.m reads a mat-file that contains:
%   x         = p*n array of input data (p data points; n variables)
%               The data in each of the 'p' rows must have been measured at 
%               the same time instant for all the 'n' variables
%   variables = n-vector of strings with the names of the variables
%
% [The name of the mat-file containing the input data is entered in section 
% PROBLEM DEFINITION. Additional user input may be entered in that section
% and in the following (OPTIONS); otherwise the algorithm is 
% executed with the default settings]                                                    
% 
% Outputs: Results are saved in a mat-file which contains the options 
% (options), the names of the variables (variables), the data points (x), 
% the name of the data file (datafile), and a structure called 'Output'
% with the following fields:
%   Output.MI         = n*n*(nlags+1) array, mutual information (several lags)
%   Output.MIl        = normalized mutual information (Linfoot 1957)
%   Output.MIm        = normalized mutual information (Michaels et al 1998)
%   Output.MIs        = normalized mutual information (Studholme et al 1999)
%   Output.H1         = n-vector of entropies
%   Output.H2         = n*n*(nlags+1) array, joint entropy
%   Output.H3         = n*n*n array of joint entropy among 3 variables
%   Output.H4         = n*n*n*n array of joint entropy among 4 variables
%   Output.MI3        = n*n*n array of three-way mutual information
%   Output.dist       = n*n*(nlags+1) array of distance between variables
%   Output.taumin     = n*n array of the time lags that minimize EMCdist
%   Output.cond_entr2 = n*n array of conditional entropies, H(Xm|Xn)
%   Output.cond_entr3 = n*n*n array of cond. entropies of 3 variables
%   Output.cond_entr4 = n*n*n*n array of cond. entropies of 4 variables
%   Output.con_array  = n*n array of connections between variables
%   Output.adaptThres = adaptive threshold value
%   Output.T          = n*n array of transfer entropies
%   Output.Y          = coordinates of the points from multidim. scaling
%
% Dependencies: This main file (runMIDER.m) calls mider.m
% -- which in turn calls estimateH2.m, estimateH3.m, and estimateH4.m -- 
% and plotResults.m; it also uses the subroutine arrow.m for visualization 
% (arrow.m was written by Erik A Johnson, see 'license.txt' for details)
%
% Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
% Created: 07/10/2011, last modified: 29/01/2015
%--------------------------------------------------------------------------
% Copyright (C) 2015  Alejandro Fernández Villaverde
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

fprintf(1,'\n \n Running MIDER...');
clear;
tStart = tic;
 

%--------------------------------------------------------------------------
% PROBLEM DEFINITION

% Matlab path of the data folder:
addpath(genpath('data'));    

% Matlab path of the results folder:
addpath(genpath('results')); 

% Input data file:
datafile = 'b2_enzyme_cat_chain'; %'b1_glycolysis';
                                  %'b2_enzyme_cat_chain';
                                  %'b3_small_chain';
                                  %'b4_irma_on_off';
                                  %'b4_irma_on_off_outliers'
                                  %'b5_mapk'; 
                                  %'b6_dream4_10_1';
                                  %'b7_dream4_100_1';
   
% Results file:
resultfile = 'results/test.mat';

load(datafile); 
npoints = size(x,1); % number of data points
ntotal  = size(x,2); % number of variables


%--------------------------------------------------------------------------
% OPTIONS 

% If you use Matlab, choose using the Statistics toolbox (=1) or not (=0):
options.useStatistics = 1;

% Data curation step:
options.correctOutliers = 1; % correct outliers (=1) or not (=0)

% Entropic parameter:
options.q = 1; % q = 1 (Boltzmann-Gibbs entropy) | q > 1 (Tsallis entropy) 

% Normalization of the mutual information (MI):
options.MItype = 'MI'; %'MI' | 'MImichaels' | 'MIlinfoot' | 'MIstudholme'  

% Adaptive estimation of MI -> Partition the joint space of X,Y so that the
% fraction of occupied bins with >= 5 points is at least = options.fraction:
options.fraction  = 0.1*(log10(npoints)-1); 
if options.fraction < 0.01, options.fraction = 0.01; end %lower bound=0.01

% Maximum time lag considered (>= 0):
options.taumax = 10;
if options.taumax > (npoints/2 -1)
    options.taumax = floor(npoints/2);
    fprintf(1,'\n Warning: the maximum time lag you have requested is too large for the number of data points.'); 
    fprintf(1,'\n => Setting the max. time lag to %d',options.taumax); 
end

% Number of entropy reduction (ERT) rounds to carry out (0, 1, 2, or 3):
options.ert_crit = 2;  

% Entropy reduction threshold. Enter a number between 0.0 and 0.2 to fix it
% manually, or enter 1 to use a value obtained from the data:
options.threshold = 1;

% Plot MI arrays (=1) or not (=0):
options.plotMI = 0;


%--------------------------------------------------------------------------
% RUN MIDER

if options.useStatistics == 1

    % Check if there are missing values in the input data:
    [row, col] = find(isnan(x));

    % If there are missing values, use TSR to complete the data set:
    if numel(row) ~= 0
        fprintf(1,'\n Warning: there are %d missing values in the input data file.',numel(row)/2);
        x = TSR(x);
        fprintf(1,'\n The missing values have been replaced with new data using the TSR method.');
    end

    % Detect possible outliers:
    [xmd,ol] = OUTLIERS(x);
    if numel(ol) ~= 0
        fprintf(1,'\n Warning: %d outlier(s) detected in the data file.',numel(ol)/2);
        fprintf(1,'\n Row(s) and column(s) of the outliers:\n');
        for i=1:numel(ol)/2
            disp(ol(i,:));
        end
        % Correct outliers, if that option has been selected:
        if options.correctOutliers == 1        
            x = TSR(xmd);
            fprintf(1,' The outliers have been replaced with new data generated with the TSR method.');
        else
            fprintf(1,' You have selected not to correct outliers.');
            fprintf(1,'\n Data points detected as outliers will be used without further processing.');
            fprintf(1,'\n If you prefer to replace them by adequate values, set options.correctOutliers = 1.');
        end
    end
    
    % Run MIDER in Matlab/Octave:
    Output = mider(x,options);

else
        
    fprintf(1,'\n You have chosen not to use functions from the Statistics toolbox.');
    fprintf(1,'\n Hence MIDER will run with limited capabilities (no data curation nor visual output).');
    fprintf(1,'\n If you have installed the Statistics toolbox and want to use it, set options.useStatistics = 1;');
    
    % Run MIDER in Matlab/Octave:
    Output = mider(x,options);

end

%--------------------------------------------------------------------------
% SAVE RESULTS

 save(resultfile,'Output','variables','x','options','datafile');


%--------------------------------------------------------------------------
% VISUALIZATION

if options.useStatistics == 1
    plotResults(Output,ntotal,variables,options)
end


%--------------------------------------------------------------------------
% TIME

fprintf(1,'\n Total wall clock time taken for MIDER: %d seconds.\n ',toc(tStart)); 
