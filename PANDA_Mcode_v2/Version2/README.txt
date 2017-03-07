README file for RunPANDA,m
Written by Kimberly Glass (kglass@jimmy.harvard.edu), available under CC BY-SA. A copy of the lisense is available at: http://creativecommons.org/licenses/by-sa/3.0/
As academic code it is provided without warranty. Please contact the author with any comments/questions/concerns.
Last updated July, 2015.

RunPANDA.m is written in the MATLAB/octave programming language and can be run by typing the name of the program at the command-prompt within MATLAB or Octave.
A the top of the script there are several parameters that can be easily set, such as the value for the update parameter (alpha) and the names for files containing the input data.

DATA:
Within the "YeastData" directory, the file "YeastNetwork_allTFxGene" contains the final z-score edge weights for the networks analyzed in "Passing Messages between Biological Networks to Refine Predicted Interactions" (Glass et. al. PLoS One, 2013). This folder also contains the expression, motif and PPI data files used to generate the networks.  Please see the supplemental material of the publication for more information on this data.

The "ToyData" directory has some small 'toy' datafiles that can be used to quickly run/understand PANDA.

FILE FORMATS:
Expression data file: In the expression data file the first column must contain gene names, and each subsequent column should contain expression values for those genes across conditions/samples (see "ToyExpressionData.txt").  Note that if fewer than three conditions are contained in the expression data file, PANDA will initialize the co-regulatory network to an identity matrix. This file can also contain header lines, so long as those lines are preceeded by a hashtag (#).

Motif data file: The motif data file contains three columns (see "ToyMotifData.txt").  The first column should contain regulators, the second those regulator's target genes, and the third, an initial weight to give the interactions (recommend 1).  Any potential regulator to gene interaction not specified in this file is given a default interaction weight of 0.  There should not be any multiply-defined edges in this file and all genes in this file must be included in the expression-data file. This file cannot contain header lines.

Protein interaction file: The protein interaction file (optional) contains three columns (see "ToyPPIData.txt").  The first two columns should contain a pair of regulators, and the third, an initial weight to give the interactions (recommend 1). This file should not contain multiply-defined edges and all TFs in this file should also be contained in the motif-data file. This file cannot contain header lines.

Outputted network file contain four columns:
TF \t Gene \t Motif-prediction \t PANDA-prediction
The values in the fourth, "PANDA-prediction", column, can loosely be interpreted as Z-scores.
