KLD2D <- function(edge, bgIndex, fgIndex, expData, method=c("integers","bandwidth")[2]){
  #' Calculates the Kullback Leibler divergence (KLD) for the 2D gene expression of X and Y in two conditions 
  #' @param edge a vector holding the names of the two genes for which KLD will be calculated
  #' @param fgIndex a vector of length M holding the indices of the samples in X and Y in condition 1
  #' @param bgIndex a vector of length M holding the indices of the samples in X and Y in condition 2
  #' @param expData a numeric matrix holding the expression data
  #' @docType methods
  #' @return The KL divergence value
  
  
  require(KernSmooth)
  x <- rank(expData[edge[1],c(bgIndex,fgIndex)], ties.method="average")
  y <- rank(expData[edge[2],c(bgIndex,fgIndex)], ties.method="average")
  N <- length(x)
  bgI <- 1:length(bgIndex)
  fgI <- length(bgIndex)+(1:length(fgIndex))
  
  # to calculate the Gaussian kernel smoothing we first need to estimate the bandwidth
  fgWidth <- c(bw.nrd(x[fgI]),bw.nrd(y[fgI]))
  bgWidth <- c(bw.nrd(x[bgI]),bw.nrd(y[bgI]))
  
  gridSize <- switch(method,
                     integers  = c(N, N),
                     bandwidth = ceiling(N / c(min(fgWidth[1], bgWidth[1]), 
                                               min(fgWidth[2], bgWidth[2]))))
  # gridSize <- pmax(gridSize,5) # make sure there are at least 25 points in total
  ranges <- list(x=c(1, N), y=c(1, N))
  
  ## debug
  #   print(cbind(x[fgI],y[fgI]))
  #   print(fgWidth)
  #   print(ranges)
  #   print(gridSize)
  
  fgSmooth <- bkde2D(x=cbind(x[fgI], y[fgI]), bandwidth=fgWidth, range.x=ranges, gridsize=gridSize)
  fgP <- fgSmooth$fhat
  bgSmooth <- bkde2D(x=cbind(x[bgI], y[bgI]), bandwidth=bgWidth, range.x=ranges, gridsize=gridSize)
  bgP <- bgSmooth$fhat
  
  # make sure there are no zeros in the smooth function (since we will take a log of that)
  #fgP[fgP==0] <- min(fgP[fgP>0])/100
  fgP <- pmax(fgP, 1e-20)
  #bgP[bgP==0] <- min(bgP[bgP>0])/100
  bgP <- pmax(bgP, 1e-20)
  
  fgP <- fgP / sum(fgP)
#   contour(fgSmooth$x1,fgSmooth$x2,fgP)
#   points(x[fgI],y[fgI],pch=16)
  
  bgP <- bgP / sum(bgP)
#   contour(bgSmooth$x1,bgSmooth$x2,bgP)
#   points(x[bgI],y[bgI],pch=16)
  return((sum(fgP * log(fgP / bgP)) + 
          sum(bgP * log(bgP / fgP))) / 2)
}