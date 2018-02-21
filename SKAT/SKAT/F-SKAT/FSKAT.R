########################################################
# SKAT for Non-continuous Traits in Familial GWAS Data #
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


FSKAT <- function(phenotype, genotypes, id, fa, mo, family="binomial", weights=NULL, covariates=NULL, acc=1e-6, append.write=NULL){
library(MASS)
library(mgcv)

# Regular checks
n<-length(phenotype)
if(!is.data.frame(genotypes)) stop("genotypes should be a data.frame!")
if(ncol(genotypes)!=n+2) stop("Number of individuals inconsistent between phenotype and genotypes. Check your data...")
if(length(id)!=n) stop("Number of individuals inconsistent between phenotype and id. Check your data...")
if(!is.null(covariates)) {
  if(nrow(covariates)!=n) stop("Number of individuals inconsistent between phenotype and covariates. Check your data...")}

y <- as.matrix(phenotype)
K <- kinship(id, fa, mo)

genotypes2 <- split(genotypes, genotypes[,1])

# Weight according to beta density function
if(is.null(weights)){
W <- list()
for (k in 1:length(genotypes2)){
 w <- dbeta(colMeans(t(genotypes2[[k]][,-c(1:2)]))/2, 1, 25)
 W[[k]] <- diag(w^2)}}
else if(!is.null(weights)){
weights2 <- split(weights, weights[,1])
W <- list()
for (k in 1:length(weights2)){
 W[[k]] <- diag(weights2[[k]][,-c(1:2)])
}}

intercept <- rep(1,length(id))
if(is.null(covariates)){
   X <- as.matrix(intercept)
   exprs<-paste("y ~ 1")}
else if(!is.null(covariates)){
   X <- cbind(intercept, as.matrix(covariates))
   X <- as.matrix(X)
   exprs<-paste("y ~ ", paste(names(covariates), collapse=" + "))}

 cs.K <- corSymm(2*K[lower.tri(K)],fixed=T)
 id <- as.factor(id)
 id <- as.matrix(id)
 cs.K <- Initialize(cs.K,data=id)
 data_pre <- data.frame(id=id,y=y)

 if(is.null(covariates)) data <- data_pre
 else if(!is.null(covariates)) data <- cbind(data_pre, covariates)

 model <- glmmPQL2(as.formula(exprs), random=~1|id, correlation=cs.K, data=data, family=family, control=lmeControl(opt="optim"))

 beta <- model$fit$coefficients$fixed
 V <- extract.lme.cov(model$fit, data)
 V_inv <- solve(V)

 # test beta value
 #beta2 <- solve(t(X)%*%V_inv%*%X)%*%t(X)%*%V_inv%*%as.matrix(model$y_star)
 
 res = as.matrix(model$y_star - X%*%beta)

Q<-eig<-evals<-tmpout<-list()
p<-vector()
 P = V - X %*% solve(t(X) %*% V_inv %*% X) %*% t(X)
for (k in 1:length(genotypes2)){  
 G = t(genotypes2[[k]][,-c(1:2)])
 Q[[k]] = t(res) %*% V_inv %*% G %*% W[[k]] %*% t(G) %*% V_inv %*% res
 eig[[k]] = eigen(sqrt(W[[k]]) %*% t(G) %*% V_inv %*% P %*% V_inv %*% G %*% sqrt(W[[k]]), symmetric=T, only.values=T)
 evals[[k]] = eig[[k]]$values[eig[[k]]$values>1e-6*eig[[k]]$values[1]]

 tmpout[[k]]<-davies(Q[[k]], evals[[k]], acc=acc)
 p[k]<-tmpout[[k]]$Qq
 if(!is.null(append.write)){
  write.table(t(c(names(genotypes2)[k], signif(p[k], digits=6))), file=append.write, row.names=F, col.names=F, append=T, quote=F) }
}
cbind(names(genotypes2), signif(p, digits=6))
}

# References:
# Davies RB. 1980. The distribution of a linear combination of chi-square random variables. Journal of the Royal Statistical Society.Series C (Applied Statistics) 29:323-333.
# Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X. 2011. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89:82-93.

