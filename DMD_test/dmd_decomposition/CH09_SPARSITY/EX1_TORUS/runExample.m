clear all, close all, clc
addpath ../utils/
% PARAMETERS
n = 128;
K = 5; % degree of sparsity
T = 2; % duration of integration
dt = 0.01; % time-step
r = 10; % number of PCA modes to keep
noisemag = 0.0; % magnitude of additive white noise
filename = 'DATA/RUN';
saveFLAG = 0;
movieFLAG = 0;
getParms % generate Damping rate, etc.

%% GENERATE SPARSE DATA
saveFLAG=0;
noisemag = 0.;
getSparseData

%% PLOT MOVIE
plotData

%% COMPUTE FFT MODES
computeFFTModes

%% COMPUTE POD MODES
computePODModes

%% COMPUTE DMD MODES
computeDMDModes

%% PROJECT DATA
projectData

%% COMPUTE DMD MODES USING COMPRESSIVE SENSING
compressedDMD