\name{printDeMAND}
\alias{printDeMAND}
\title{Basic methods for class demand}
\usage{
printDeMAND(x)
}
\arguments{
\item{x}{An instance of class demand which includes: a gene expression data, annotation information, and a molecular interaction network.}
}
\value{printDeMAND returnssummary information about the diggit object}
\description{This document lists a series of basic methods for the class DeMAND}
\examples{
data(inputExample)
dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
printDeMAND(dobj)
}



