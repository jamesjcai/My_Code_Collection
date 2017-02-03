library(glmnet)
cisgenos <- readRDS(file="exampledata_predictable/cisgenos.rds")
exppheno <- readRDS(file="exampledata_predictable/exppheno.rds")
groupid <- readRDS(file="exampledata_predictable/groupid.rds")

n_k_folds<-10;
alpha<-1;

fit <- cv.glmnet(as.matrix(cisgenos),
                 as.vector(exppheno),
                 nfolds = n_k_folds,
                 alpha = alpha,
                 keep = TRUE,
                 foldid = groupid,
                 parallel = FALSE)


# Pull info from fit to find the best lambda   
fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))
# Needs to be min or max depending on cv measure (MSE min, AUC max, ...)
best.lam <- fit.df[which.min(fit.df[,1]),]
cvm.best <- best.lam[,1]
lambda.best <- best.lam[,2]
# Position of best lambda in cv.glmnet output
nrow.best <- best.lam[,3]
# Get the betas from the best lambda value
ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
ret[ret == 0.0] <- NA
# Pull the non-zero betas from model
as.vector(ret[which(!is.na(ret)),])


pred.mat <- fit$fit.preval[,nrow.best]
res <- summary(lm(exppheno~pred.mat))
rsq <- res$r.squared
pval <- res$coef[2,4]





# ------------

bestbetas <- tryCatch(
  { 
    fit2 <- cv.glmnet(as.matrix(cisgenos),
                     as.vector(exppheno),
                     nfolds = n_k_folds,
                     alpha = alpha,
                     keep = TRUE,
                     foldid = groupid,
                     parallel = FALSE)
    # Pull info from fit to find the best lambda   
    fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))
    # Needs to be min or max depending on cv measure (MSE min, AUC max, ...)
    best.lam <- fit.df[which.min(fit.df[,1]),]
    cvm.best <- best.lam[,1]
    lambda.best <- best.lam[,2]
    # Position of best lambda in cv.glmnet output
    nrow.best <- best.lam[,3]
    # Get the betas from the best lambda value
    ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
    ret[ret == 0.0] <- NA
    # Pull the non-zero betas from model
    as.vector(ret[which(!is.na(ret)),])
    
    
  },
  error = function(cond) {
    # Should fire only when all predictors have 0 variance.
    message('Error with gene ' %&% gene %&% ', index ' %&% i)
    message(geterrmessage())
    return(data.frame())
  }
)



