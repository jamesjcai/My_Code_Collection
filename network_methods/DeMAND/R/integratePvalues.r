integratePvalues <- function(g, network, expData, ppi, keepLeaves){
  #' This function integrates the p-values of all the edges surrounding a gene using Fisher's method,
  #' and uses Brown's method to correct for correlations between the p-values.
  #' @param g the name of the gene
  #' @param network a matrix or data frame holding the gene names in the first two columns, followed by the KLD value and followed by the p-balue
  #' @param expData gene expression data with gene names as rownames
  #' @docType methods
  #' @return the corrected p-value for gene g

  # find all the interactions involving g
  neighborEdges<- which(network[,1] %in% g | network[,2] %in% g)
  N <- length(neighborEdges)
  if(N<2)
    if (keepLeaves) {
      return(as.numeric(network[neighborEdges,4][1]))
    } else {
      return(1)
    }
  
  # calculate the Fisher method integrated chi square value
  pvals <- as.numeric(as.vector(network[neighborEdges,4]))
  pvals <- pmax(pvals,1e-20)
  fisherChisq <- -2*sum(log(pvals))
  
  # get the expression of the neighbors to calculate correlations
  neighborGenes <- apply(network[neighborEdges,1:2],1,function(x) setdiff(x,g))
  if(length(unique(neighborGenes))<2)
    if (keepLeaves) {
      return(as.numeric(network[neighborEdges,4][1]))
    } else {
      return(1)
    }
  neighborExp <- expData[neighborGenes,]
  gExp <- expData[g,]
    
  # since we are trying to correlate esges, we estimate the residuals from a fit to gene g
  resids <- lm(t(neighborExp)~gExp)$residuals
  # however, for protein-protein interactions we use the original data
  resids[ ,ppi[neighborEdges]] <- t(expData[neighborGenes, ][ppi[neighborEdges],])
  covMat <- cov(resids)
  
  # calculate the correction to Fisher's method
  covMat <- covMat[lower.tri(covMat)]
  I <- covMat<0 # Brown's method requires a different correction for positive and negative values
  covMat[I] <- covMat[I]*(3.27 + 0.71*covMat[I])
  covMat[!I] <- covMat[!I]*(3.25 + 0.75*covMat[!I])
  # the correction coefficient is
  c <- 4*N/(4*N + 2*sum(covMat))
  dFreedom <- 2*N*c
  correctedChisq <- fisherChisq*c
  
  return(pchisq(q=correctedChisq,df=dFreedom,lower.tail=F))
}
