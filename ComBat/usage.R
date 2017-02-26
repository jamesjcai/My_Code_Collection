library(R.matlab)
source('ComBat.R')

# skip=1 if the first column contains probeid
outdata<-ComBat('infile.txt','infile2.txt',skip=0,write=F)

writeMat('outfile.mat', data=as.matrix(outdata))