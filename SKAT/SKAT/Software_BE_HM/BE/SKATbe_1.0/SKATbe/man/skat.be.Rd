\name{skat.be}
\alias{skat.be}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SKAT test with backward elimination and subsampling.
}
\description{
Either SKAT or a burden test is used to test for association between a set of genetic variants and a phenotypic state. Backwards elimination is used to find an optimal set of associated variants. Subsampling is used to speed computation.
}
\usage{
skat.be( Z, y, w, File.Out, ro=0, N.SIMR=300, basedir=getwd() )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
	A matrix of genotype data. Each row corresponds to an individual and each column corresponds to a variant. Variant information is coded as 0 (homozygous major allele), 1 (heterozygous), or 2 (homozyous minor allele).
}
  \item{y}{
	A vector of phenotype values.
}
  \item{w}{
	A vector of weights to be applied to the variants. The length of the vector should be equal to the number of columns in Z.
}
  \item{File.Out}{
	Name of directory to save result in. The directory will be created if it does not already exist.
}
  \item{ro}{
	Set the type of test to be used. Set to 0 for SKAT and 1 for burden. Default is 0.
}
  \item{N.SIMR}{
  	The number of subsamples to use. Default is 300.
}
  \item{basedir}{
  	The base directory for the file arguments. Default is the current working directory.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	Produces a list of candidate variants and saves them to files specified by File.Out. The files are:
	
	File.Out_BIG.ro1 - Gives the full results from the backwards elimination, with the variant eliminated in each step in the first column and the resulting pvalue in the second column.
	
	File.Out_BIG.ro1.curr_set.ro1 - The final set of variants selected in the backwards elimination.
	
	File.Out.ro1.big_curr_set.ro1 - The final set of variants selected in the backwards elimination for each random subsample of the data. One subsample per line.

	File.Out.ro1.big_set20.ro1 - The set of variants selected in each subsample. One subsample per line.	
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
# Load example data
data(VPS13B)

skat.be(z.vps13b,y.vps13b,w.vps13b,File.Out="ExampleResults")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SKAT }
\keyword{ burden }% __ONLY ONE__ keyword per line
