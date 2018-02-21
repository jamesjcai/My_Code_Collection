\name{skat.cl}
\alias{skat.cl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SKAT test with backward elimination and subsampling.
}
\description{
Either SKAT or a burden test is used to test for association between a set of genetic variants and a phenotypic state. Backwards elimination is used to find an optimal set of associated variants. Subsampling is used to speed computation.
}
\usage{
skat.cl( File.Ped, File.W, File.Out, ro=0, N.SIMR=300, basedir=getwd() )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{File.Ped}{
	File name specifying a PED file with genotype and phenotype information. Each row corresponds to an individual. The first six columns contain family id, individual id, paternal id, maternal id, sex, and phenotype, and each remaining column corresponds to a variant. Variant information is coded as 0 (homozygous major allele), 1 (heterozygous), or 2 (homozyous minor allele). The functions in SKATbe only use the phenotype and genotype data, but assumes the other columns are present. If data is missing for columns 1-5, it may be filled in with a place holder character (e.g., "-" or ".").
}
  \item{File.W}{
	File name containing a vector of weights to apply to variants. Weights should be non-negative numbers, one per line in the file. The number of weights should match the number of genotype columns in the PED file.
}
  \item{File.Out}{
	The file to write results to. Several files are produced, with suffixes added to File.Out for the complete file name.
}
  \item{ro}{
	Set the type of test to be used. Set to 0 for SKAT and 1 for burden. Default is 0.
}
  \item{N.SIMR}{
  	The number of subsamples to use. Default is 300.
}
  \item{basedir}{
  	The base directory for the file arguments. Default is the current working directory. If basedir does not exist it will be created.
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

# From the command line run the following in the directory containing RunSKATbe.R
# R --slave --args VPS13B.ped weights_VPS13B_nonsyn.txt VPS13B_results ExampleData 1 2000 < RunSKATbe.R

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SKAT }
\keyword{ burden }% __ONLY ONE__ keyword per line
