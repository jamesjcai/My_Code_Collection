#'
#'  This function is based on the realization that drugs affect the protein activity of their targets, but not necessarily their mRNA expression levels. In contrast, the change in protein activity directly affects the mRNA expression levels of downstream genes. Based on this hypothesis, DeMAND identifies drug MoA by comparing gene expression profiles following drug perturbation with control samples, and computing the change in the individual interactions within a pre-determined integrated transcriptional and post-translational regulatory model (interactome). For each edge in the interactome we determine the two-dimensional probability distribution of the gene expression levels both in the control state, and following drug treatment. Any changes in the probability distribution are estimated using the Kullback-Leibler (KL) divergence, from which we determine the statistical significance of the dysregulation of each edge. In the second step of DeMAND, we interrogate each gene independently to determine whether its interactions are enriched in dysregulated ones, suggesting that it is a candidate mechanism of action.
#'
#' @param dobj Instance of class demand
#' @param fgIndex Sample indices of Drug treated samples
#' @param bgIndex Sample indices of DMSO treated samples
#' @param verbose Whether to print progress text
#' @param method How to choose the bins to evaluate the KL-divergence. The options are using intervals based on the kernel bandwidth (default) or on integer points in the data
#' @param keepLeaves Whether to provide a p-value for genes the have only a single neighbor in the network (default is FALSE)
#' @param alpha The cutoff for estimating a p-value using pareto fitting (default=0.05)
#' @return Objet of class demand with updated moa slot
#' @docType methods
#' @examples
#' data(inputExample, package="demand")
#' dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
#' dobj <- runDeMAND(dobj, fgIndex=caseIndex, bgIndex=controlIndex)
#' show(dobj)
#' @export


runDeMAND <- function (x, fgIndex=NULL, bgIndex=NULL, verbose=TRUE, 
                       method="bandwidth", keepLeaves=FALSE, alpha=0.05){
  
  # Parameter check
  if (is.null(fgIndex) | is.null(bgIndex)) 
    stop("Please provide sample (column of expression data) indices of case/control samples")
  if (any(is.na(c(fgIndex, bgIndex)))){
    warning('Case indices contain NA values. These values are ignored')
    fgIndex <- as.vector(na.exclude(fgIndex))
    bgIndex <- as.vector(na.exclude(bgIndex))
  }
  if (length(fgIndex) < 3 | length(bgIndex) < 3) 
    stop("The number of samples in each class should be at least three")
  if (length(fgIndex) < 6 | length(bgIndex) < 6) 
    warning("DeMAND requires six samples (in each class) for optimal performance")
  expData <- x@exp
  if (any(is.na(expData))) 
    stop("Expression data contains NA values")
  if (any(is.infinite(x@exp))) 
    warning("Expression data contains infinite values")
  annot <- x@anno[ ,2]
  if (any(is.na(annot))) 
    warning("Annotation data contains NA values")
  inputNetwork <- x@network
  if (any(is.na(inputNetwork[, 1:2]))){
    warning("The network contains NA values, removing those lines")
    inputNetwork <- inputNetwork[apply(inputNetwork, 1, function(x) !any(is.na(x))), ]
    x@network <- inputNetwork
  } 
  platformGene <- unique(annot)
  
  vmsg <- function(x, verb=verbose){
    if (verb)
      message(x)
  } 
  
  # To reduce the networks using genes appear in the expression data
  vmsg("Pruning the network")
  edgesToKeep <- apply(inputNetwork[, 1:2], 1, function(gg) all(gg %in% platformGene))
  interactome <- inputNetwork[edgesToKeep, ]
  rm(edgesToKeep)
  analGene <- intersect(unique(as.vector(interactome[, 1:2])), platformGene)
  # remove duplicate lines (for instance x-y and y-x)
  interactome[, 1:2] <- t(apply(interactome[,1:2], 1, function(gg) c(min(gg),max(gg))))
  dups <- duplicated(interactome)
  interactome <- interactome[!dups, 1:2]
  ppi <- if ("ppi" %in% colnames(inputNetwork)) {
    as.numeric(inputNetwork[!dups, "ppi"])==1 
  } else {
    rep(F, sum(!dups))
  }
  
  # To preprocess the expression data. if there are multiple probes for one gene, 
  # select one with the maximum coefficient of variation (CV)
  vmsg("Keeping best probe per gene")
  CV <- sqrt(rowMeans(expData ^ 2) - rowMeans(expData) ^ 2) / rowMeans(expData)
  oGenes <- order(annot, -CV) # the minus sign is to sort from high to low in the CV
  dups <- duplicated(annot[oGenes])
  expData <- expData[oGenes[!dups], ]
  row.names(expData) <- annot[oGenes[!dups]]
  expData <- expData[analGene, ]
  
  # load the function that calculates 2D KL-divergence
  #source('KLD2D.r')  
  
  # get statistical significance of the KLD value using null distribution from the randomized network
  getKLDpvalue <- function(kld, nullKLD) {
    rs <- sum(nullKLD >= kld) / length(nullKLD)
    return(min(1, rs))
  }
  
  vmsg("Make a null distribution for KL divergence.....")
  # create a random network of the same size as the original network
  # (but no less than 1K and no more than 10K)
  p1 <- sample(analGene, min(max(length(analGene), 1e3), 1e4), replace = T)
  p2 <- sample(analGene, min(max(length(analGene), 1e3), 1e4), replace = T)
  pKeep <- !(p1 == p2)
  permuteInteractome <- cbind(p1[pKeep], p2[pKeep])
  permuteInteractome <- t(apply(permuteInteractome, 1, function(x) c(min(x), max(x))))
  dups <- duplicated(permuteInteractome)
  permuteInteractome <- permuteInteractome[!dups, ]
  
  # get null distribution of KLD values
  nullBgIndex <- bgIndex#sample(x=c(bgIndex,fgIndex),size=length(bgIndex),replace=F)
  nullFgIndex <- fgIndex#setdiff(c(bgIndex,fgIndex),nullBgIndex)
  nullKLD <- apply(permuteInteractome, 1, KLD2D, nullBgIndex, nullFgIndex, expData, method)
  
  vmsg("Measure dysregulation of the interactions.....")
  
  # get KLD for all edges
  KLDmat <- apply(interactome, 1, KLD2D, bgIndex, fgIndex, expData, method)
  
  # get pvalue (of KLD) for all the edges using pareto distribution of the tails
  # source('pareto.R')
  pfit <- pareto.fit(data=nullKLD, threshold=quantile(nullKLD, probs=1 - alpha))
  KLDpvec <- ppareto(x=KLDmat, threshold=pfit$xmin, exponent=pfit$exponent, lower.tail=F) * alpha
  # use the counts of the nullKLD to estimate higher p-values
  pToReplace <- KLDmat < quantile(nullKLD, probs=1 - alpha)
  KLDpvec[pToReplace] <- sapply(X=KLDmat[pToReplace], FUN=function(K) getKLDpvalue(K, nullKLD))
  edgeKLD <- cbind(interactome, KLDmat, KLDpvec)
  colnames(edgeKLD) <- c("gene1", "gene2", "KLD", "KLD.p")
  
  # To combine p-values
  vmsg("Estimate dysregulation of the genes.....")
  # source('integratePvalues.r')
  intPval <- apply(as.matrix(analGene), 1, integratePvalues, edgeKLD, expData, ppi, keepLeaves)
  
  intPvalAdjustedp <- p.adjust(intPval, "fdr")
  finalPrediction <- data.frame(moaGene = analGene, Pvalue = intPval, FDR = intPvalAdjustedp)
  finalPrediction <- finalPrediction[order(intPval, decreasing = F), ]
  rownames(finalPrediction) <- 1:dim(finalPrediction)[1]
  
  # To update moa slot in the demand object
  x@moa <- finalPrediction
  x@KLD <- as.data.frame(edgeKLD[order(KLDpvec), ])
  return(x)
}

