\docType{class}
\name{demand-instance}
\alias{demand-instance}
\title{The DeMAND instance}
\description{
This instance stores parameters and results of theDeMAND algorithm
}
\arguments{
\item{exp}{A N by M numeric matrix and the rows are N probes and the columns are M samples.}
\item{anno}{A N by 2 character matrix. The rows are probes but the order should be the same with the \emph{demand} matrix.  The first column can be anything (usually probe IDs) but the second column should
 includes Official Gene Symbol information for each probe.}
\item{network}{A K by 4 character matrix which contains K interactions. The 1st column and the 2nd column contain pairs of interacted genes. The 3rd and 4th columns indicate whether the interactions are pr
otein-protein interaction (ppi) or protein-DNA (pdi) interaction. Column name should be as follows: c("Gene1", "Gene2", "ppi", "pdi")}
\item{moa}{A data frame contains DeMAND MoA predictions (e.g. Gene, p-value, adjusted p-values)}
\item{KLD}{A matrix containing the KL-divergence of the interactions that were analysed, the KL-divergence that was evaluated, and the p-value associated with the divergenece}
}

