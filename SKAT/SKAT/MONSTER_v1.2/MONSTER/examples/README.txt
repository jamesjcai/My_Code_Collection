# This directory provides example input files.

# To test if the pre-compiled binary executable works on your system:
../bin/MONSTER -p pheno_ex.txt -g geno_ex.txt -s SNP_ex.txt -k kin_ex.txt
# To test the binary executable compiled by yourself, simply replace ../bin/ by the path under which the executable file is created

# This command should produce two output files: MONSTER.out and MONSTER.param
# MONSTER.out should look like
SNP_set_ID	n_individual	n_SNP rho_MONSTER	p_MONSTER
SNPset1	769	3	0.0	0.050095966928
SNPset2	766	2	0.0	0.865208281368

# MONSTER.param should look like
SNPset1 
Estimate	nullMLE	SE_nullMLE
Heritability	0.449534	0.064429
Additive_Var	5.421304	0.949546
Error_Var	6.638524	0.688633
Intercept	6.728266	0.409844
Covariate_1	0.557100	0.231901

SNPset2 
Estimate	nullMLE	SE_nullMLE
Heritability	0.466789	0.064509
Additive_Var	5.636473	0.962558
Error_Var	6.438509	0.692624
Intercept	6.734958	0.409984
Covariate_1	0.558878	0.231408

# More example commands and detailed explanations can be found in the documentation.

