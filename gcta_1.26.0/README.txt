PROGRAM: GCTA

DESCRIPTION: Genome-wide Complex Trait Analysis

AUTHOR: Jian Yang, Hong Lee, Mike Goddard and Peter Visscher

CONTACT: jian.yang@uq.edu.au

YEAR: 2010-2013

MIT License

DOCUMENTATION: http://www.complextraitgenomics.com/software/gcta/

INSTALLATION: When you have download a zip or gzipped archive with an
executable binary, no installation is necessary.

USAGE: Type "gcta64" or "./gcta64" from the command line followed by the
options of choice (see documentation) NOTE: you probably need to run 
"chmod a+x gcta" to get the correct permission to execute the program.

EXAMPLE DATA: Four example files test.bed, test.bim, test.fam and test.phen 
are included in the distribution; for example, once GCTA is installed try running:

     gcta64 --bfile test --make-grm --out test

     gcta64 --reml --grm test --pheno test.phen --out test

     etc...
