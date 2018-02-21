\name{scatit}
\alias{scatit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SKAT test with backward elimination.
}
\description{
Either SKAT or a burden test is used to test for association between a set of genetic variants and a phenotypic state. Backwards elimination is used to find an optimal set of associated variants.
}
\usage{
scatit(Z, y, File.Out, ro)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
	A genotype matrix. Each row corresponds to an individual, and each column corresponds to a variant. Variant information is coded as 0 (homozygous major allele), 1 (heterozygous), 2 (homozyous minor allele), or 9 (missing).
}
  \item{y}{
	A vector specifying the phenotype of each individual.
}
  \item{File.Out}{
	The file to write results to.
}
  \item{ro}{
	Set the type of test to be used. Set to 0 for SKAT and 1 for burden.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Produces a list of candidate variants and saves them to File.Out	
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
