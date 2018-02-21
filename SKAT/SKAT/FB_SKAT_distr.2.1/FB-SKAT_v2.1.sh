#!/bin/bash

#!/bin/bash

###############################################################################
#
# Below are the main input parameters. You need to
# Change the values according to your own dataset.
# No change of the CPP code is required

# In our environment it compiles with the command:
# g++ $(gsl-config --cflags) -o FB-SKAT_v2.1 FB-SKAT_v2.1.cpp $(gsl-config --libs)

###############################################################################



pedigreeFile=data.ped #read pedigree file; 9 is missing data; 1 is affected; 0 unaffected; data.ped
variantPassFile=variant_pass.txt #pass information for each variant; if this file is not available, set all to 1.; variant_pass.txt
genesFile=genes.txt  #read gene file : gene names & start/end positions; genes.txt
weightsFile=weights.txt #read weights for variants; e.g. 1 for nonsynonymous and 0 otherwise; if this information is missing, set all weights to 1. weights.txt
resultFile=results_ # results files
mendelianErrorsFile=mendelian_errors_  #mendelian errors files
npermut=10000 #number of random permutations
MAF=0.05 #maximum MAF on the variants to be tested
mend_th=0.01 #maximum mendelian error rate at a marker that is allowed;
ro=0 #specifies which test to run: ro=0 SKAT; ro=1 BURDEN.
min_mark=3 #minimum number of markers in a gene to be analyzed.

## Run SKAT test (ro=0)
./FB-SKAT_v2 ${pedigreeFile} ${variantPassFile} ${genesFile} ${weightsFile} ${resultFile} ${mendelianErrorsFile} ${npermut} ${MAF} ${mend_th} 0 ${min_mark}
## Change ro to "1" to run BURDEN test
./FB-SKAT_v2 ${pedigreeFile} ${variantPassFile} ${genesFile} ${weightsFile} ${resultFile} ${mendelianErrorsFile} ${npermut} ${MAF} ${mend_th} 1 ${min_mark}

echo "DONE"
