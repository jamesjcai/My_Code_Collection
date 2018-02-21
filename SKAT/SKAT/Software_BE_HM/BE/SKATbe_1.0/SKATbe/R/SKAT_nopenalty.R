#!/bin/sh

#SKAT(Z, obj, kernel = "linear.weighted",
#       method="davies", weights.beta=c(1,25), weights=NULL,
#       impute.method="fixed", r.corr=0, is_check_genotype=TRUE,
#       is_dosage = FALSE, missing_cutoff=0.15 )


#Z: a numeric genotype matrix with each row as a different
#   individual and each column as a separate gene/snp. Each
#   genotype should be coded as 0, 1, 2, and 9 (or NA) for AA,
#   Aa, aa, and missing, where A is a major allele and a is a
#   minor allele. Missing genotypes will be imputed by the simple
#   Hardy-Weinberg equilibrium (HWE) based imputation.
 # obj: an output object of the SKAT_Null_Model function.

###
#  To Run
#  R --slave -f SKAT_iuliana.R "--args cc100.ped SKAT1 0 1 10"
###

#rho=0 - SKAT
#rho=1 - Burden
#pD.SKAT<-SKAT(Z, objD, r.corr=0)$p.value
#pD.Burden<-SKAT(Z, objD, r.corr=1)$p.value


library(SKAT)

scatit<-function(Z, y, File.Out, ro){

objD<-SKAT_Null_Model(y ~ 1, out_type="D", Adjustment=FALSE)

out.p<-NULL


flag=1
x<-dim(Z)[2]
curr_set=1:x;

while (flag==1) {

result <- try({

pD<-SKAT(Z, objD, r.corr=ro)$p.value
porig=pD;

i<-0
out.p<-rbind(out.p,c(i, pD))

x<-dim(Z)[2]
cat ("Z", dim(Z), "\n" )

if (x==2) {flag=0; break;} 
pd1=vector('numeric',length=x);

for (i in 1:x){
    pD<-SKAT(Z[,-i], objD, r.corr=ro)$p.value
    pd1[i]=pD;
    out.p<-rbind(out.p,c(i,pD));
}

#calculate smallest p value by removing one i
minp=min(pd1);

if (minp>porig) { flag=0; }
px=sort(pd1,index.return=TRUE);
i0=px$ix[1];
Z=as.matrix(Z[,-i0]);
if (flag==1) curr_set=curr_set[-i0];

colnames(out.p)<-c("IndRemoved", "pD")

write.table(out.p, file = File.Out, append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = ".", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"))


})  #ERROR HANDLING
if(class(result) == "try-error") next;

}## Close while loop

logfile<-paste(File.Out, '.curr_set', '.ro', ro, sep='')
write(t(curr_set),file=logfile,ncolumns=1);
#write(curr_set,file=logfile,ncolumns=length(curr_set));
}

