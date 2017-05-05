%--------------------------------------------------------------------------
% runMIDER: Call the MIDER algorithm (mider.m), save results and plot them
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
% Dependencies: This main file (runMIDER.m) calls mider.m -- which in turn 
% calls estimateH2.m, estimateH3.m, and estimateH4.m -- and plotResults.m
% runMIDER.m also uses the subroutine arrow.m for visualization 
% (arrow.m was written by Erik A Johnson, see 'license.txt' for details)
%
% Written by Alejandro Fernández Villaverde (afvillaverde@iim.csic.es)
% Created: 07/10/2011, last modified: 17/01/2014
%--------------------------------------------------------------------------
% Copyright (C) 2013  Alejandro Fernández Villaverde
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

clear;
fprintf(1,'\n \n Running MIDER...'); 


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
                                  %'b5_mapk'; 
                                  %'b6_dream4_10_1';
                                  %'b7_dream4_100_1';
   
% Results file:
resultfile = 'results/results_mider_b2.mat';

load(datafile); 
npoints = size(x,1); % number of data points
ntotal  = size(x,2); % number of variables


%--------------------------------------------------------------------------
% OPTIONS 

% Entropic parameter:
options.q = 1; % q = 1 (Boltzmann-Gibbs entropy) | q > 1 (Tsallis entropy) 

% Normalization of the mutual information (MI)
options.MItype = 'MI'; %'MI' | 'MImichaels' | 'MIlinfoot' | 'MIstudholme'  

% Adaptive estimation of MI -> Partition the joint space of X,Y so that the
% fraction of occupied bins with >= 5 points is at least = options.fraction:
options.fraction  = 0.1*(log10(npoints)-1); 
if options.fraction < 0.01, options.fraction = 0.01; end %lower bound=0.01

% Maximum time lag considered (> 0):
options.taumax = 10;

% Number of entropy reduction (ERT) rounds to carry out (0, 1, 2, or 3):
options.ert_crit = 2;  

% Entropy reduction threshold. Enter a number between 0.0 and 0.2 to fix it
% manually, or choose 'adapt' to use a value obtained from the data:
options.threshold = 'adapt';

% Plot MI arrays (=1) or not (=0):
options.plotMI = 1;


%--------------------------------------------------------------------------
% RUN MIDER

Output = mider(x,options);


%--------------------------------------------------------------------------
% SAVE RESULTS

save(resultfile,'Output','variables','x','options','datafile');


%--------------------------------------------------------------------------
% VISUALIZATION

plotResults(Output,ntotal,variables,options)
