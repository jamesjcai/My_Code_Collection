# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scater", version = "3.8")
# BiocManager::install("scran", version = "3.8")
library(scran)
library(scater)

X<-read.csv('X.txt',header = FALSE)
Y<-read.csv('Y.txt',header = FALSE)
Z1<-read.csv('Z1.txt',header = FALSE)
Z2<-read.csv('Z2.txt',header = FALSE)

makeSCE <- function(X){
  mName <- deparse(substitute(X))
  X <- as.matrix(X)
  rownames(X) <- paste0("G",seq_len(nrow(X)))
  colnames(X) <- paste0(mName,seq_len(ncol(X)))
  X <- SingleCellExperiment(list(counts=X))
  X <- computeSumFactors(X)
  X <- normalize(X)
  return(X)
}

X <- makeSCE(X)
Y <- makeSCE(Y)
Z1 <- makeSCE(Z1)
Z2 <- makeSCE(Z2)

# https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#7_batch_correction

out <- fastMNN(X, Y)
out
dim(out$corrected)
