\name{demandClass}
\alias{demandClass}
\title{The demand class constructor}
\usage{
demandClass(exp, anno, network, moa=NULL, KLD=NULL)
}
\arguments{
\item{exp}{A N-by-M numeric matrix where the rows represent N probes (or genes) and the columns represent M samples.}
\item{anno}{A N-by-2 character matrix where the rows represent probes or genes in the same order as the \code{exp} matrix. The first column must hold the probe id or gene name as appears in the \code{exp} matrix, and the second column should hold their corresponding names (e.g gene symbol) as appears in the \code{network} matrix}
\item{network}{A K-by-L (L>1) character matrix containing K interactions. The 1st column and the 2nd column contain the names of the interacting genes. If the following columns include a column called "ppi" then genes connected by ppi will be used without estimating the residuals in the Brown correction method. The rest of the columns are used only to distinguish duplicates.}
\item{moa}{Filled my the \code{runDeMAND} function. A data frame containing DeMAND MoA predictions (i.e. Gene, p-value, adjusted p-values)}
\item{KLD}{Filled my the \code{runDeMAND} function. A matrix containing the interactions that were analysed, their KL-divergence, and the p-value associated with the divergenece}
}
\value{
Instance of class demand
}
\description{
This function generates demand class instances
}
\examples{
## Load toy example
data(inputExample)
dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
printDeMAND(dobj)
}

