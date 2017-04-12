library(SKAT)
library(R.matlab)
D<-readMat("skatinfile.mat");
attach(D);
# source('skatinfile.r')
obj1<-SKAT_Null_Model(expr ~ 1, out_type="C")
pskat1<-SKAT(g012, obj1, method="optimal.adj")$p.value;


obj2<-SKAT_Null_Model(expr.residual ~ 1, out_type="C")
pskat2<-SKAT(g012, obj2, method="optimal.adj")$p.value;

writeMat("skatoutfile.mat",pskat1=pskat1,pskat2=pskat2)
