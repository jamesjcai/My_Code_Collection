MAIN FILES

computePOD.m  — Computes POD for cylinder wake data
computeDMD.m  — Computes DMD for cylinder wake data

Both files load data using:

    load ../DATA/FLUIDS/CYLINDER_ALL.mat

This data could also be generated using the IBPM code and loaded using the following support files.  

SUPPORT FILES

loadIBPM.m — loads velocity field snapshot from IBPM code
loadDATA.m — loads multiple velocity field snapshots and creates “VORTALL” data file, saved in ../DATA/FLUIDS/CYLINDER_ALL.mat

CCcool.mat — Color scheme to plot vorticity
