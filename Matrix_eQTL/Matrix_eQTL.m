% Matrix eQTL by Andrey A. Shabalin
% http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
%

clear;clc;
% Load the library
addpath Matrix_eQTL_Matlab

%% Settings

% Linear model to use, modelANOVA or modelLINEAR
useModel = modelLINEAR; % modelANOVA or modelLINEAR

% Genotype file name
SNP_file_name = 'Sample_Data/SNP.txt';

% Gene expression file name
expression_file_name = 'Sample_Data/GE.txt';

% Covariates file name
% Set to [] for no covariates
covariates_file_name = [];  %'Sample_Data/Covariates.txt';

% Output file name
output_file_name = 'Sample_Data/eQTL_results_M.txt';

% Only associations significant at this level will be saved
pvOutputThreshold = 1e-2;

% Error covariance matrix
% Set to [] for identity.
errorCovariance = [];


%% Load genotype data

snps = SlicedData;
snps.fileDelimiter = char(9); % the TAB character
snps.fileOmitCharacters = 'NA'; % denote missing values;
snps.fileSkipRows = 1; % one row of column labels
snps.fileSkipColumns = 1; % one column of row labels
snps.fileSliceSize = 2000; % read file in pieces of 2,000 rows
snps.LoadFile( SNP_file_name );

%% Load gene expression data

gene = SlicedData;
gene.fileDelimiter = char(9); % the TAB character
gene.fileOmitCharacters = 'NA'; % denote missing values;
gene.fileSkipRows = 1; % one row of column labels
gene.fileSkipColumns = 1; % one column of row labels
gene.fileSliceSize = 2000; % read file in pieces of 2,000 rows
gene.LoadFile( expression_file_name );

%% Load covariates

cvrt = SlicedData;
cvrt.fileDelimiter = char(9); % the TAB character
cvrt.fileOmitCharacters = 'NA'; % denote missing values;
cvrt.fileSkipRows = 1; % one row of column labels
cvrt.fileSkipColumns = 1; % one column of row labels
cvrt.fileSliceSize = 2000; % read file in one piece
if(~isempty(covariates_file_name))
	cvrt.LoadFile( covariates_file_name );
end;

%% Run the analysis

Matrix_eQTL_engine( snps, ...
                    gene, ...
                    cvrt, ...
                    output_file_name, ...
                    pvOutputThreshold, ...
                    useModel, ...
                    errorCovariance, ...
                    true);
