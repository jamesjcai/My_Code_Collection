library(R.matlab)
library(SKAT)
data(SKAT.example)
#names(SKAT.example)
attach(SKAT.example)

obj<-SKAT_Null_Model(y.b ~ 1, out_type="D")
#SKAT(Z, obj)$p.value

obj<-SKAT_Null_Model(y.c ~ X, out_type="C")
SKAT(Z, obj)$p.value


library(AssotesteR)


SKATx <-
  function(y, X, kernel="linear", weights=NULL, a=1, b=25, perm=NULL)
  {
    ## checking argumetns
    if (!is.vector(y) || mode(y) != "numeric")
      stop("argument 'y' must be a numeric vector")
    if (any(is.na(y))) 
      stop("Sorry =(   No missing data allowed in argument 'y' ")	
    if (!all(y %in% c(0, 1)))
      stop("Sorry =(   argument 'y' must contain only 0 and 1")
    if(!is.matrix(X) & !is.data.frame(X))
      stop("argument 'X' must be a matrix or data.frame")
    if (nrow(X) != length(y)) 
      stop("'X' and 'y' have different lengths")
    if (!is.matrix(X)) X = as.matrix(X)
    if (!(kernel %in% c('linear', 'wlinear', 'quadratic', 'IBS', 'wIBS', 'twowayx')))
      stop(paste("\n", "Sorry =(   I don't recognize kernel:", kernel, "\n",
                 "Choose one from: 'linear', 'wlinear', 'quadratic', 'IBS', 'wIBS', 'twowayx'")) 
    if (!is.null(weights))
    {
      if(length(weights) != ncol(X)) 
        stop("length of 'weights' does not match number of columns in 'X'")
      if(!is.vector(weights) || mode(weights) != "numeric" || any(weights) <= 0)
        stop("argument 'weights' must be a numeric vector with positive values")
    } 
    if (mode(a) != "numeric" || length(a) != 1 || a <= 0)
      stop("argument 'a' must be a positive number")
    if (mode(b) != "numeric" || length(b) != 1 || b <= 0)
      stop("argument 'b' must be a positive number")
    if (!is.null(perm))
    {
      if (mode(perm) != "numeric" || length(perm) != 1 || perm < 0 || (perm %% 1) !=0) 
      {
        warning("Argument 'perm' incorrectly defined. Value perm=NULL is used")
        perm = NULL
      }
    } else perm=0
    
    ## how many obs and variants
    n = nrow(X)
    p = ncol(X)
    
    ## get weights for linear-weighted or IBS-weighted
    if (kernel=="wlinear" || kernel=="wIBS" && is.null(weights))
    {   
      # MAF across cases and controls combined
      MAF <- colMeans(X, na.rm=TRUE) / 2
      # weights of the variants
      weights = dbeta(MAF, a, b)
    }	
    
    ## apply kernel
    K = switch(kernel,
               "linear" = X %*% t(X),
               "wlinear" = ((X %*% diag(weights^2)) %*% t(X) ), 
               "quadratic" = (X %*% t(X) + 1) ^ 2,
               "IBS" = kernel_IBS(X, n, p),
               "wIBS" = kernel_wIBS(X, n, p, weights),
               "twowayx" = kernel_twowayx(X, n, p)
    )
    print(K)
    ## score statistic Q (follows a Chi-square distr)
    y.new = y - mean(y)
    Q = t(y.new) %*% K %*% y.new / 2 
    
    ## parameters to estimate distr of Q    
    mu = rep(mean(y), n)
    V = diag(mu * (1-mu))
    ones = rep(1, n)
    P = V - V %*% ones %*% solve(t(ones) %*% V %*% ones) %*% t(ones) %*% V
    PK = P %*% K
    muQ = sum(diag(PK)) / 2
    
   
    itau = sum(PK * t(PK)) / 2
    
    itausig = sum(PK * t(P)) / 2
    
    
    
    isigma = sum(P^2) / 2
    
    iest = itau - ((itausig^2)/isigma)
    skale = iest / (2 * muQ)
    print(skale)
    
    
    degfr = 2 * (muQ^2) / iest
    print(skale)
    print(degfr)
    ## p-value
    skat.stat = as.numeric(Q)
    asym.pval = 1 - pchisq(skat.stat/skale, df=degfr)
    
    ## permutations
    perm.pval = NA
    if (perm > 0)
    {
      x.perm = rep(0, perm)
      for (i in 1:perm)
      {
        perm.sample = sample(1:length(y))
        y.perm = y[perm.sample]
        y.new = y.perm - mean(y.perm)
        Q.perm = t(y.new) %*% K %*% y.new / 2
        x.perm[i] = as.numeric(Q.perm)
      }
      # p-value 
      perm.pval = sum(x.perm > skat.stat) / perm
    } else perm = "NULL"   
    
    ## results
    name = "SKAT: Sequence Kernel Association Test"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), perm, kernel)
    names(arg.spec) = c("cases", "controls", "variants", "n.perms", "kernel")
    res = list(skat.stat = skat.stat, 
               asym.pval = asym.pval, 
               perm.pval = perm.pval, 
               args = arg.spec, 
               name = name)
    class(res) = "assoctest"
    return(res)
  }



SKATx(y.b,Z,kernel="IBS")
