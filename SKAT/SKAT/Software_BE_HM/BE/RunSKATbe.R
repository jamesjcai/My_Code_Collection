##################
#
# To run execute the command
# R --slave --args <file.ped> <file.weights> <dir.save> <base> <ro> <N.SIMR> < RunSKATbe.R
#
# Example:
# R --slave --args AAK1.ped weights_AAK1.txt AAK1_results ExampleData 1 300 < RunSKATbe.R
#
# Note: AAK1_results is a directory and not a file. It will be created if it
# does not already exist.
#
# See the documentation for skat.cl() in the SKATbe package for more options.
#
##################


library(SKATbe)

args <- commandArgs(trailingOnly=TRUE)

ped <- args[1]
weights <- args[2]
save <- args[3]
base <- args[4]
r <- args[5]
N <- args[6]

skat.cl( ped, weights, save, basedir=base, ro=r, N.SIMR=N )


