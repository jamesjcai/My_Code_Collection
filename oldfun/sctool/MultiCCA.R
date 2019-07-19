library(Seurat)
library(Matrix)

# Preparing inputs
X<-read.csv('X.txt',header = FALSE)
Y<-read.csv('Y.txt',header = FALSE)
Z1<-read.csv('Z1.txt',header = FALSE)
Z2<-read.csv('Z2.txt',header = FALSE)
prepareData <- function(X){
  mName <- deparse(substitute(X))
  X <- as.matrix(X)
  rownames(X) <- paste0("G",seq_len(nrow(X)))
  colnames(X) <- paste0(mName,seq_len(ncol(X)))
  X <- CreateSeuratObject(X)
  X <- NormalizeData(X)
  X <- FindVariableGenes(X, do.plot = F, display.progress = F)
  X <- Seurat::ScaleData(X)
  X@meta.data$tech <- mName
  return(X)
}

X <- prepareData(X)
Y <- prepareData(Y)
Z1 <- prepareData(Z1)
Z2 <- prepareData(Z2)

# Generating the list
ob.list <- list(X, Y, Z1, Z2)


# Run multi-set CCA
pancreas.integrated <- RunMultiCCA(ob.list, num.ccs = 15)

# CC Selection
MetageneBicorPlot(pancreas.integrated, grouping.var = "tech", dims.eval = 1:15)

# Run rare non-overlapping filtering
pancreas.integrated <- CalcVarExpRatio(object = pancreas.integrated, reduction.type = "pca",
                                       grouping.var = "tech", dims.use = 1:10)
pancreas.integrated <- SubsetData(pancreas.integrated, subset.name = "var.ratio.pca",
                                  accept.low = 0.5)

# Alignment
pancreas.integrated <- AlignSubspace(pancreas.integrated,
                                     reduction.type = "cca",
                                     grouping.var = "tech",
                                     dims.align = 1:10)

# t-SNE and Clustering
pancreas.integrated <- FindClusters(pancreas.integrated, reduction.type = "cca.aligned",
                                    dims.use = 1:10, save.SNN = T, resolution = 0.4)
pancreas.integrated <- RunTSNE(pancreas.integrated,
                               reduction.use = "cca.aligned",
                               dims.use = 1:10)

# Visualization
TSNEPlot(pancreas.integrated, do.label = T)

