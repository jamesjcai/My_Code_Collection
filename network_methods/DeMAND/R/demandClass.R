#'The DeMAND class
#'
#'This class stores parameters and results of the DeMAND algorithm
#'  
#'@section Slots:
#'    \describe{
#'      \item{\code{exp}:}{Matrix containing the gene expression data}
#'      \item{\code{anno}:}{Matrix containing annotation for the gene expression data}
#'      \item{\code{network}:}{Matrix containing molecular interaction network}
#'      \item{\code{moa}:}{Data frame containing a list of genes and their DeMAND p-values}
#'      \item{\code{KLD}:}{Matrix containting the molecular interactions that were used for the analysis plus the KL-divergence and the p-value that was calculated}
#'}
#'@rdname demand-instance
setClass("demand", representation(exp="matrix", anno="matrix", network="matrix", moa="data.frame", KLD="data.frame"))

#' The demand class constructor
#' 
#' This function generates demand class objects
#' 
#' @param exp Numeric matrix with features in rows and samples in columns
#' @param anno Character matrix containing annotation for the gene expression data
#' @param network Character matrix containing molecular interaction network
#' @param moa Data.frame containing a list of genes and their DeMAND p-values
#' @param KLD a matrix contating molecular interactions (gene names), their KL-divergence, and the DeMAND p-value
#' @return Object of class demand
#' @examples
#' data(inputExample)
#' dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
#' @export
demandClass <- function(exp, anno, network, moa=NULL, KLD=NULL) {
    if (is.null(exp)) stop("This class requires an expression data (numeric matrix) to be generated.")
    if (is.null(anno)) stop("This class requires an annotation of the expression data (character matrix) to be generated.")
    if (is.null(network)) stop("This class requires a network (character matrix) to be generated.")
    if (is.null(moa)) moa=data.frame()
    if (is.null(KLD)) KLD=data.frame()
    new("demand", exp=as.matrix(exp), anno=as.matrix(anno), network=as.matrix(network), moa=moa, KLD=KLD)
}
