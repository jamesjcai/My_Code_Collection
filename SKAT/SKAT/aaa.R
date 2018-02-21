library(SKAT)
data(SKAT.example)
names(SKAT.example)
attach(SKAT.example)

obj<-SKAT_Null_Model(y.b ~ X, out_type="D")

SKAT(Z, obj)$p.value

# library(AssotesteR)


my_check <-
  function(y, X, perm)
  {
    ## Internal function to check argument:s y, X, perm
    
    # y as numeric vector
    if (!is.vector(y) || mode(y) != "numeric")
      stop("Sorry, argument 'y' must be a numeric vector")
    # no misssing values in y
    if (any(is.na(y))) 
      stop("No missing data allowed in argument 'y' ")	
    # binary values (0, 1) in y
    if (!all(y %in% c(0, 1)))
      stop("Argument 'y' must contain only 0 and 1")
    # X as matrix or data frame
    if(!is.matrix(X) & !is.data.frame(X))
      stop("Argument 'X' must be a matrix or data.frame")    
    # compatibility between X and y
    if (nrow(X) != length(y)) 
      stop("'X' and 'y' have different lengths")
    # force X as matrix
    if (!is.matrix(X)) X = as.matrix(X)
    # permutations
    if (mode(perm) != "numeric" || length(perm) != 1
        || perm < 0 || (perm %% 1) !=0) 
    {
      warning("argument 'perm' incorrectly defined. Value perm=100 is used")
      perm = 100
    }
    # results
    list(y=y, X=X, perm=perm)
  }



my_cmc_method <-
  function(casecon, X.new)
  {
    # Internal function for CMAT method
    ## number of individuals N, cases nA, controls nU
    N = nrow(X.new)
    nA = sum(casecon)
    nU = N - nA
    ## matrix of genotypes in cases
    Xx = X.new[casecon==1,]  
    ## matrix of genotypes in controls  
    Yy = X.new[casecon==0,] 
    ## get means
    Xx.mean = colMeans(Xx, na.rm=TRUE)
    Yy.mean = colMeans(Yy, na.rm=TRUE)
    ## center matrices Xx and Yy
    Dx = sweep(Xx, 2, Xx.mean)
    Dx
    Dy = sweep(Yy, 2, Yy.mean)
    
    ## pooled covariance matrix
    if (sum(complete.cases(X.new)) == N)  # no missing values
    {  
      COV = (t(Dx) %*% Dx + t(Dy) %*% Dy) / (N-2)
    } else {  # with missing values
      ## covariance matrix of cases
      tDx = t(Dx)
      Sx = matrix(0, ncol(X.new), ncol(X.new))
      for (i in 1:nrow(tDx))
      {
        for (j in i:ncol(Dx))
        {
          Sx[i,j] = sum(tDx[i,] * Dx[,j], na.rm=TRUE)
        }
      }
      sx.diag = diag(Sx)
      Sx = Sx + t(Sx)
      diag(Sx) = sx.diag
      ## covariance matrix of controls
      tDy = t(Dy)
      Sy = matrix(0, ncol(X.new), ncol(X.new))
      for (i in 1:nrow(tDy))
      {
        for (j in i:ncol(Dy))
        {
          Sy[i,j] = sum(tDy[i,] * Dy[,j], na.rm=TRUE)
        }
      }
      sy.diag = diag(Sy)
      Sy = Sy + t(Sy)
      diag(Sy) = sy.diag
      ## pooled covariance matrix
      COV = (1/(N-2)) * (Sx + Sy)	
    }
    
    ## general inverse
    if (nrow(COV) == 1) # only one variant
    { 
      if (COV < 1e-8) COV = 1e-8
      COV.inv = 1 / COV
    } else {
      COV.eigen = eigen(COV)
      eig.vals = COV.eigen$values  
      inv.vals = ifelse(abs(eig.vals) <= 1e-8, 0, 1/eig.vals)
      EV = solve(COV.eigen$vectors)
      COV.inv = t(EV) %*% diag(inv.vals) %*% EV
    }	
    
    ## Hotellings T2 statistic
    stat = t(Xx.mean - Yy.mean) %*% COV.inv %*% (Xx.mean - Yy.mean) * nA * nU / N
    as.numeric(stat)
  }



CMC <-
  function(y, X, maf=0.05, perm=100)
  {
    ## checking arguments
    Xy_perm = my_check(y, X, perm)
    y = Xy_perm$y
    X = Xy_perm$X
    perm = Xy_perm$perm
    
    ## number of individuals N
    N = nrow(X)
    ## get minor allele frequencies
    MAF = colMeans(X, na.rm=TRUE) / 2   
    ## how many variants < maf
    rare.maf = MAF < maf
    rare = sum(rare.maf)
    ## collapsing
    if (rare <= 1) 
    {   
      # if rare variants <= 1, then NO collapse is needed
      X.new = X
    } else {
      # collapsing rare variants into one column
      X.collaps = rowSums(X[,rare.maf], na.rm=TRUE)
      X.collaps[X.collaps != 0] = 1
      # joining collapsed to common variants
      X.new = cbind(X[,!rare.maf], X.collaps)	   
    }
    ## change values to -1, 0, 1
    X.new = X.new - 1
    ## number of new variants
    M = ncol(X.new)
    ## Hotellings T2 statistic
    cmc.stat = my_cmc_method(y, X.new)
    
    ## Asymptotic p-values
    # under the null hypothesis T2 follows an F distribution 
    f.stat = cmc.stat * (N-M-1)/(M*(N-2))
    df1 = M          # degrees of freedom  
    df2 = N - M - 1  # degrees of freedom  
    asym.pval = 1 - pf(f.stat, df1, df2)
    
    ## under the alternative hyposthesis T2 follows a chi-square distr
    # pval = 1 - pchisq(cmc.stat, df=M)
    
    ## permutations
    perm.pval = NA
    if (perm > 0)
    {
      x.perm = rep(0, perm)
      for (i in 1:perm)
      {
        perm.sample = sample(1:length(y))
        x.perm[i] = my_cmc_method(y[perm.sample], X.new) 
      }
      # p-value 
      perm.pval = sum(x.perm > cmc.stat) / perm
    }
    
    ## results
    name = "CMC: Combined Multivariate and Collapsing Method"
    arg.spec = c(sum(y), length(y)-sum(y), ncol(X), rare,  maf, perm)
    arg.spec = as.character(arg.spec)
    names(arg.spec) = c("cases", "controls", "variants", "rarevar", "maf", "perm")	
    res = list(cmc.stat = cmc.stat, 
               asym.pval = asym.pval, 
               perm.pval = perm.pval, 
               args = arg.spec, 
               name = name)
    class(res) = "assoctest"
    return(res)
  }

CMC(y.b,Z)

