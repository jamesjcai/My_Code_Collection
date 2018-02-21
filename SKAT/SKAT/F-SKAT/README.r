########################################################
# SKAT for Binary Traits in Familial GWAS Data         #
# Qi Yan (qi.yan@chp.edu)                              #
# Version 1.0                                          #
########################################################

# NOTES:
# This function (F-SKAT) is used to perform SKAT analysis (Wu et al., 2011) for non-continuous traits in familial GWAS data
# In the final correlation matrix, the covariance between parent and offspring is 0.5, the covariance between siblings is also 0.5
# It requires R package "kinship" to get parent-offspring kinship matrix
# After downloading the package, you can install it using
# R CMD INSTALL -l ~/R_library_directory kinship_1.1.3.tar.gz
# library(kinship, lib.loc="~/R_library_directory")

# It requires the R package "CompQuadForm" to compute the p-value using Davies' method (Davies, 1980)
# After downloading the package, you can install it using
# R CMD INSTALL -l ~/R_library_directory CompQuadForm_1.4.tar.gz
# library(CompQuadForm, lib.loc="~/R_library_directory")

# It also requires "glmmPQL.s" that is modified glmmPQL function used to fit null model
# source("/working directory/glmmPQL.s")

# FUNCTION ARGUMENTS:
# phenotype: A vector of quantitative trait in the analysis (class: vector). The order should match the vector id. No missing.
# genotypes: 1st column: gene name; 2nd column: snp name; 3rd-end columns: A matrix of genotypes for each subject (class: data.frame). The order of 3rd-end columns should match id. Coded as 0, 1, 2 and no missing.
# id: A vector of id (class: vector). It can be either numeric or character. The id indicates each subject. Make sure it is not factor. No missing.
# fa: A vector of father id (class: vector). It can be either numeric or character. The father id indicates the father of each subject. If this subject has no father in this data, the value is set to "NA". Make sure it is not factor.
# mo: A vector of mother id (class: vector). It can be either numeric or character. The mother id indicates the mother of each subject. If this subject has no mother in this data, the value is set to "NA". Make sure it is not factor.
# family: Type of phenotype. (Default="binomial")
# covariates: A matrix of covariates (class: data.frame). The order of rows should match the vector id. Default NULL. No missing.
# weights: 1st column: gene name; 2nd column: snp name; 3rd column: A vector with the length equal to the number of variants in the test (class: data.frame). Default is Null indicating weight function as described in Wu's SKAT paper
# acc: Accuracy of numerical integration used in Davies' method. Default 1e-6.
# append.write: The name of pvalue output file. Write out p-values in real time. Don't need to wait until all genes are processed to get the final output.

# OUTPUT:
# pvalue: non-continuous trait family SKAT (F-SKAT) p-value


### For testing multiple genes (GWAS)
library(kinship, lib.loc="~/R_library_directory")
library(CompQuadForm, lib.loc="~/R_library_directory")
setwd('/working directory')
source("glmmPQL.s")
source("FSKAT.R")

# Subject IDs are numeric
y <- read.table(file="phenotype.txt", header=T)
gene <- read.table(file="gene.txt", header=T)
weights <- read.table(file="weights.txt", header=T)
covariates <- read.table(file="covariates.txt", header=T)
pvalue1 <- FSKAT(phenotype=y$y, genotypes=gene, id=y$id, fa=y$fa, mo=y$mo, family="binomial", covariates=NULL, weights=NULL, append.write="./pvalues.out")

pvalue1 <- FSKAT(phenotype=y$y, genotypes=gene, id=y$id, fa=y$fa, mo=y$mo, family="binomial", covariates=NULL, weights=NULL)
pvalue2 <- FSKAT(phenotype=y$y, genotypes=gene, id=y$id, fa=y$fa, mo=y$mo, family="binomial", covariates=NULL, weights=weights)
pvalue3 <- FSKAT(phenotype=y$y, genotypes=gene, id=y$id, fa=y$fa, mo=y$mo, family="binomial", covariates=covariates, weights=weights)

# Subject IDs are character
y <- read.table(file="phenotypeX.txt", header=T)
gene <- read.table(file="geneX.txt", header=T)
pvalue1 <- FSKAT(phenotype=y$y, genotypes=gene, id=as.character(y$id), fa=as.character(y$fa), mo=as.character(y$mo), family="binomial", covariates=NULL, weights=NULL)
