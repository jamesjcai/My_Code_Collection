\name{runDeMAND}
\alias{runDeMAND}
\title{Run DeMAND algorithm}
\usage{
runDeMAND(x, fgIndex, bgIndex, verbose=TRUE, method="bandwidth", 
          keepLeaves=FALSE, alpha=0.05)
}
\arguments{
\item{x}{An instance of class demand which includes: a gene expression data, annotation information, and a molecular interaction network.}
\item{fgIndex}{A numeric vector contains indices of columns which represent case samples (e.g. drug treated). The sample size should be greater than 3.}
\item{bgIndex}{A numeric vector contains indices of columns which represent control samples (e.g. drug treated). The sample size should be greater than 3.}
\item{verbose}{A boolean value (TRUE by defalut) indicating whether to print progression outputs}
\item{method}{A string value indicating whether to evaluate the KL-divergence using grid points based on \code{'bandwidth'} (default) or on the \code{'integers'} space}
\item{keepLeaves}{A boolean value indicating whether to return a p-value for genes the have only a single neighbor in the network (default is FALSE, which returns 1 for such genes)}
\item{alpha}{The cutoff for estimating a p-value using pareto fitting (default=0.05)}
}
\value{
Returns a DeMAND class object holding the same exp, anno, and network slots as the input, and where the moa and KLD slots hold the results of DeMAND algorithm. The moa slot is a matrix containing a list of genes, corresponding p-value, and adjusted p-value. The KLD slot is a matrix with the gene names of the edges that were evaluated, their KL-divergence, and the p-value assigned to that divergence.
}
\description{
DeMAND.

This function is based on the realization that drugs affect the protein activity of their targets, but not necessarily their mRNA expression levels. In contrast, the change in protein activity directly affects the mRNA expression levels of downstream genes. Based on this hypothesis, DeMAND identifies drug MoA by comparing gene expression profiles following drug perturbation with control samples, and computing the change in the individual interactions within a pre-determined integrated transcriptional and post-translational regulatory model (interactome).
}
\details{
For each edge in the interactome we determine the two-dimensional probability distribution of the gene expression levels both in the control state, and following drug treatment. Any changes in the probability distribution are estimated using the Kullback-Leibler (KL) divergence, from which we determine the statistical significance of the dysregulation of each edge. In the second step of DeMAND, we interrogate each gene independently to determine whether its interactions are enriched in dysregulated ones, suggesting that it is a candidate mechanism of action.
}
\examples{
## Load toy example
data(inputExample)
dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork, moa=NULL)
dobj <- runDeMAND(dobj, fgIndex=caseIndex, bgIndex=controlIndex)
## results (head)
printDeMAND(dobj)
## results (all)
print(dobj@moa)
print(dobj@KLD)
}
\keyword{Action}
\keyword{Drug}
\keyword{Mechanism}
\keyword{of}

