#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <Eigen>
#include "MatVec.h"
#include "MONSTER_read.h"
#include "MONSTER_null.h"
#include "MONSTER_test.h"


using namespace Eigen;
typedef Map<Matrix<double, Dynamic, Dynamic, RowMajor> > MapRowXd;
typedef Map<Matrix<double, Dynamic, Dynamic> > MapColXd;
typedef Map<Matrix<int, Dynamic, Dynamic, RowMajor> > MapRowXi;
typedef Map<Matrix<int, Dynamic, Dynamic> > MapColXi;

#define MAXLEN 30
#define MAXSNPSET 1000


int findSNP(int n_allsnp, char **allsnp, char snp[]);
void get_eigen(int n_fam, struct FAMILY *family, double **U, double *lambda);
void get_eigen2(int n_fam, struct FAMILY *family, int *n_indiv, double ***U, double **lambda);
void imputation(int n_fam, int n_snp, struct FAMILY *family, int BLUP, int BLUE, int AVE);
void formatdata(int n_fam, int n_cov, int n_snp, int n_indiv, struct FAMILY *family, double **y, double ***x, double ***G);
void formatdata_null(int n_fam, int n_cov, int n_indiv, struct FAMILY *family, double **y, double ***x);





int main(int argc, char **argv){

    printf("\n"
           "--------------------------------------------------------------------------------\n"
           "|                                  MONSTER                                     |\n"
		   "| (SNP-set Association Testing for Quantitative Traits in Related Individuals) |\n"
           "|                         Version 1.2 - June 26, 2015                          |\n"
           "|                                                                              |\n"
           "|                          Copyright(C) 2014, 2015                             |\n"
           "|                      Duo Jiang and Mary Sara McPeek                          |\n"
           "|                                                                              |\n"
           "|                                 Homepage:                                    |\n"
           "|             http://galton.uchicago.edu/~mcpeek/software/MONSTER              |\n"
           "--------------------------------------------------------------------------------\n\n"
           );

	int n_grid = 11;
	char pedfile[MAXLEN] = "pheno.txt", genofile[MAXLEN]="geno.txt", snpfile[MAXLEN]="SNP.txt", kinfile[MAXLEN]="kin.txt", eigfile[MAXLEN], mcut_str[MAXLEN];
	char outfile[]="MONSTER.out", paramfile[]="MONSTER.param", residfile[]="MONSTER.resid", eigoutfile[]="MONSTER.eig";
	int arg, nm=0, pfile=0, gfile=0, sfile=0, kfile=0, erecyc=0, efile=0, comparison=0, mc=0, BLUP=1, BLUE=0,res=0, AVE=0, w;
	double mcut = 0, **U, *lambda, *weight, weight_tmp, *y, **x, **G, pvalue, rho_optim, pvalue_burden, pvalue_second;
	int n_indiv_all, n_cov, n_fam, i, n_indiv, n_allsnp, *snpstart, 
		tmp, length=MAXSNPSET*MAXSNPLEN+MAXSNPLEN+3, n, n_snp, *snpsetstart, skip, maxfamsize, j;
	struct FAMILY *family, *family_tmp;
	char **allsnp, snpset[MAXSNPLEN], line[MAXSNPSET*MAXSNPLEN+MAXSNPLEN+3], line2[MAXSNPSET*MAXSNPLEN+MAXSNPLEN+3], snp_tmp[MAXSNPLEN], snpset2[MAXSNPLEN];
	FILE *snpf, *outf, *paramf, *residf, *eigoutf;
	struct NullResult *nullResult;

	if (argc > 1) {
		for (arg = 1; arg < argc && argv[arg][0] == '-'; arg++) {
			switch (argv[arg][1]) {
			case 'n':
				nm = 1;
				break;
			case 'p':
				strncpy(pedfile, argv[++arg], MAXLEN);
				pfile = 1;
				break;
			case 'g':
				strncpy(genofile, argv[++arg], MAXLEN);
				gfile = 1;
				break;
			case 's':
				strncpy(snpfile, argv[++arg], MAXLEN);
				sfile = 1;
				break;	
			case 'k':
				strncpy(kinfile, argv[++arg], MAXLEN);
				kfile = 1;
				break;	
			case 'c':
				comparison = 1;
				break;
			case 'r':
				res = 1;
				break;
			case 'm':
				strncpy(mcut_str, argv[++arg], MAXLEN);
				mcut = atof(mcut_str);				
				mc = 1;
				break;	
			case 'E':
				BLUE = 1;
				BLUP = 0;
				break;		
			case 'A':
				AVE = 1;
				BLUP = 0;
				break;					
			case 'e':
				strncpy(eigfile, argv[++arg], MAXLEN);
				erecyc = 1;
				efile = 1;
				break;					
			default:
				printf("Unknown option \"%s\"\n", argv[arg]);
				exit(1);
			}
		}
	}

	if (!pfile) {
		printf("Pedigree and phenotype file: %s\n", pedfile);
	} else{
		printf("User specified pedigree and phenotype file: %s\n", pedfile);
	}
	if (!kfile) {
		printf("Kinship file: %s\n", kinfile);
	} else{
		printf("User specified kinship file: %s\n", kinfile);
	}
	if(nm){
		printf("The program will only fit the null model, without performing any association tests.\n");
		if(res==1){
			printf("Phenotypic residuals requested under the null model.\n");
		}
		if(gfile){
			printf("WARNING: Genotype file %s is ignored!\n", genofile);
		}
		if(sfile){
			printf("WARNING: SNP list file %s is ignored!\n", snpfile);		
		}
		if(comparison==1){
			printf("WARNING: The flag -c is ignored because of the option -n.\n");
		}
		if(mc==1 && mcut!=1){
			printf("WARNING: The flag -m is ignored because of the option -n.\n");
		}
		mcut=1;
		if(BLUE==1){
			printf("WARNING: The flag -E is ignored because of the flag -n!\n");
		}
		if(AVE==1){
			printf("WARNING: The flag -A is ignored because of the flag -n!\n");
		}
	} else {
		if (!gfile) {
			printf("Genotype file: %s\n", genofile);
		}else{
			printf("User specified genotype file: %s\n", genofile);
		}
		if (!sfile) {
			printf("SNP list file: %s\n", snpfile);
		} else{
			printf("User specified SNP list file: %s\n", snpfile);
		}
		if(res==1){
			printf("Phenotypic residuals requested under the null models.\n");
		}
		if(comparison==1){
			printf("famSKAT and famBT will be performed, in addition to the MONSTER testing method.\n");
		}
		if(!mc){
			printf("Individuals with missing genotypes are excluded from the analysis.\n");
		} else{
			if(mcut<0 || mcut>1){
				printf("ERROR: Missingness threshold must be between 0 and 1!\n");
				exit(1);
			} 
			if(mcut<1){
				printf("Individuals with over %lf%% of the genotypes missing are excluded from the analysis. \n", mcut*100);
			} else {
				printf("All individuals are included in he analysis, regardless of missingness.\n");
			}
		}
		if(BLUE==1){
			if(AVE==1){
				printf("ERROR: Options -B and -A cannot be used simultaneously!\n");
				exit(1);
			} else if(mcut>0){
				printf("BLUE imputation will be used for missing genotypes, if any.\n");
			}
		} else if(AVE==1){
			if(mcut>0) printf("Sample average imputation will be used for missing genotypes, if any.\n");
		} else{
			if(mcut>0) printf("BLUP imputation will be used for missing genotypes, if any.\n");
		}
		if(erecyc==1 && mcut<1){
			printf("WARNING: the user-specified eigen-decomposition file %s is ignored, because the mssingness threshold is below 1.\n", eigfile);
			erecyc=0;
		}
	}


	count_ped(pedfile, &n_indiv_all, &n_cov, &n_fam);
    printf("\nWARNING: User should make sure these counts are consistent with the pedigree and phenotype file. \n");
    printf("n_family = %d\tn_individual = %d\tn_covariate = %d\n\n",n_fam, n_indiv_all, n_cov);
	
	family = (struct FAMILY*) malloc(n_fam * sizeof(struct FAMILY));	
	maxfamsize = read_ped(pedfile, n_indiv_all, n_fam, n_cov, family);
	if(BLUP==1 && maxfamsize>200 &&!nm) {
		printf("WARNING: With a maximum family size of %d, BLUP imputation might be slow. User may consider BLUE or sample average for faster computation.\n", maxfamsize);
	}
	read_kin (kinfile, n_fam, family);
	
	if(erecyc == 1){
		U = d2array(n_indiv_all, n_indiv_all);
		lambda = d1array(n_indiv_all);
		efile = read_eig(eigfile, n_indiv_all, U, lambda);
		n_indiv = n_indiv_all;
		if(efile == 1) {
			printf("User specified eigen-decomposition file: %s\n", eigfile);
		} else{
			printf("WARNING: Failure to find user specified eigen-decomposition file %s! The flag -e is ignored.\n", eigfile);
		}
	}
	
	if(efile == 0 && mcut == 1){
		erecyc = 1;
		U = d2array(n_indiv_all, n_indiv_all);
		lambda = d1array(n_indiv_all);
		get_eigen(n_fam, family, U, lambda);
		n_indiv = n_indiv_all;
		eigoutf = fopen(eigoutfile, "w");
		for(i=0; i<n_indiv; i++) fwrite(&(lambda[i]), sizeof(double), 1, eigoutf);
		for(i=0; i<n_indiv; i++){
			for(j=0; j<n_indiv; j++){
				fwrite(&(U[i][j]), sizeof(double), 1, eigoutf);
			}
		}
		fclose(eigoutf);
	}
	
	if(nm==1){	
		if(n_indiv <= n_cov+3){
			printf("ERROR: The null model cannot be fitted, because there are only %d individuals in the dataset!", n_indiv);
			for(i=0; i<n_fam; i++){
				family_tmp = family+i;
				free(family_tmp->indiv);
				free(family_tmp->missingrate);
				free(family_tmp->pheno);
				free_d2array(family_tmp->cov, family_tmp->famsize, n_cov+1);
				free_d2array(family_tmp->kin, family_tmp->famsize, family_tmp->famsize);
			}
			free(family);
			free(lambda);
			free_d2array(U, n_indiv, n_indiv);
			exit(1);
		}
	
		formatdata_null(n_fam, n_cov, n_indiv, family, &y, &x);
		nullResult = (struct NullResult *) malloc(sizeof(struct NullResult));
		SKAT_null (n_indiv, n_cov+1, res, y, x, lambda, U, nullResult);
		paramf = fopen(paramfile, "w");
		fprintf(paramf, "NullModel \nEstimate\tnullMLE\tSE_nullMLE\n");
		fprintf(paramf, "Heritability\t%lf\t%lf\n", (nullResult->error)[0][0], (nullResult->error)[0][1]);
		fprintf(paramf, "Additive_Var\t%lf\t%lf\n", (nullResult->error)[1][0], (nullResult->error)[1][1]);
		fprintf(paramf, "Error_Var\t%lf\t%lf\n", (nullResult->error)[2][0], (nullResult->error)[2][1]);
		fprintf(paramf, "Intercept\t%lf\t%lf\n", (nullResult->cov)[0][0], (nullResult->cov)[0][1]);
		for(i=1; i<n_cov+1; i++){
			fprintf(paramf, "Covariate_%d\t%lf\t%lf\n", i, (nullResult->cov)[i][0], (nullResult->cov)[i][1]);
		}
		fprintf(paramf, "\n");
		fclose(paramf);
		
		if(res==1){
			residf = fopen(residfile, "w");
			fprintf(residf, "NullModel \n");
			for(i=0; i<n_fam; i++){
				family_tmp = family+i;
				for(j=0; j<(family_tmp->famsize); j++){
					fprintf(residf, "%ld ", (family_tmp->indiv)[j]);
				}
			}
			fprintf(residf, "\n");
			for(i=0; i<n_indiv; i++){
				fprintf(residf, "%lf ", (nullResult->resid0)[i]);
			}
			fprintf(residf, "\n");
			for(i=0; i<n_indiv; i++){
				fprintf(residf, "%lf ", (nullResult->resid_std)[i]);
			}
			fprintf(residf, "\n");
			fclose(residf);
		}		
		for(i=0; i<n_fam; i++){
			family_tmp = family+i;
			free(family_tmp->indiv);
			free(family_tmp->missingrate);
			free(family_tmp->pheno);
			free_d2array(family_tmp->cov, family_tmp->famsize, n_cov+1);
			free_d2array(family_tmp->kin, family_tmp->famsize, family_tmp->famsize);
		}
		free(family);
		free(lambda);
		free_d2array(U, n_indiv, n_indiv);
		free(y);
		free_d2array(x, n_indiv, n_cov+1);
		
		
		free(nullResult->resid);
		free_d2array(nullResult->Sigma_inv_half, n_indiv, n_indiv);
		free_d2array(nullResult->cov, n_cov+1, 2);
		if(res==1){
			free(nullResult->resid0);
			free(nullResult->resid_std);
		}
		free(nullResult);
		printf("\nProgram completed!\n");
		if(res==1 && efile ==1){
			printf("Output files: MONSTER.param\tMONSTER.resid\n");
		} else if(res==1){
			printf("Output files: MONSTER.param\tMONSTER.resid\tMONSTER.eig\n");
		} else if(efile==1){
			printf("Output files: MONSTER.param\n");
		} else{
			printf("Output files: MONSTER.param\tMONSTER.eig\n");
		}
		return 0;
	}
	
	n_allsnp = count_allsnp(genofile);
	snpstart = i1array(n_allsnp);
	allsnp = (char**)malloc(n_allsnp * sizeof(char*));
	for (i=0; i<n_allsnp; i++) allsnp[i] = (char*)malloc(MAXSNPLEN * sizeof(char));
	read_allsnp(genofile, n_fam, family, allsnp, snpstart);
	
	if( (snpf = fopen(snpfile, "r"))==NULL) {
		printf("ERROR: failure to open the snp list file!\n");
		exit(1);
	}	
	outf = fopen(outfile, "w");
	paramf = fopen(paramfile, "w");
	if(res==1) residf = fopen(residfile, "w");
	fprintf(outf, "SNP_set_ID\tn_individual\tn_SNP rho_MONSTER\tp_MONSTER");
	if(comparison==1) fprintf(outf, "\tp_famBT\tp_famSKAT");
	fprintf(outf, "\n");
	while (!feof(snpf)){
		if(fscanf(snpf, "%s", snpset)!=1) break;
		n_snp = 0;
		skip = 0;
		fgets(line, length, snpf);
		n = sscanf(line,"%d %[^\n]", &w, line);
		if (n<2) {	
			printf("WARNING: No SNP found in SNP set %s\n", snpset);
			continue;
		}
		else {
			strcpy(line2, line);
			while (n==2){
				n=sscanf(line,"%s %[^\n]", snp_tmp, line);
				if (n>=1)(n_snp)++;		
			}
		}	
		snpsetstart = i1array(n_snp);
		for(i=0; i<n_snp; i++){
			sscanf(line2, "%s %[^\n]", snp_tmp, line2);
			snpsetstart[i] = findSNP(n_allsnp, allsnp, snp_tmp);			
			if (snpsetstart[i] <0){
				skip = 1;
				printf("WARNING: SNP %s in SNP set %s is not present in the genotype data file! This SNP set is skipped.\n", snp_tmp, snpset);
				break;
			}
			snpsetstart[i] = snpstart[snpsetstart[i]];
		}
		weight = d1array(n_snp);
		if(w==1){
			fscanf(snpf, "%s %d", snpset2, &tmp);
			if (strcmp(snpset, snpset2)!=0) {
				printf("ERROR: In SNP list input file, weights for SNPs in SNP set %s cannot be found! \n", snpset);
				exit(1);
			}
			for(i=0; i<n_snp; i++) {
				fscanf(snpf, "%lf", &weight_tmp);
				if(weight_tmp<=0){
					printf("ERROR: SNP weights must be postive! Exception found for SNP %s in SNP set %s!\n", snp_tmp, snpset);
					exit(1);
				}
				weight[i] = weight_tmp;
			}
		} else{
			for(i=0; i<n_snp; i++) weight[i] = 1;
		}
			
		if(skip==1) {
			fprintf(outf, "%s This SNP set is not tested because of an unidentified SNP.\n\n", snpset);
			fprintf(paramf, "%s \n This SNP set is not tested because of an unidentified SNP.\n\n", snpset);
			if(res==1) fprintf(residf, "%s \n This SNP set is not tested because of an unidentified SNP.\n\n", snpset);
			free(snpsetstart);
			free(weight);
			continue;
		}
	
		for(i=0; i<n_fam; i++){
			family_tmp = family+i;
			family_tmp->G_all = d2array(family_tmp->famsize, n_snp);
		}

		read_geno(genofile, n_fam, n_snp, n_indiv_all, mcut, snpsetstart, family);
		if(mcut>0) imputation(n_fam, n_snp, family, BLUP, BLUE, AVE);

		if (erecyc==0) {
			get_eigen2(n_fam, family, &n_indiv, &U, &lambda);
		}
		if(n_indiv <= n_cov+3){
			printf("WARNING: SNP set %s is not analyzed, because only %d individuals are retained based on the missingness cutoff. Unable to fit the null model!\n", snpset, n_indiv);
			fprintf(outf, "%s\t%d\t%d\tNA\tNA ", snpset, n_indiv, n_snp);
			if(comparison==1) fprintf(outf, "\tNA\tNA");
			fprintf(outf, "\n");
			fprintf(paramf, "%s \n This SNP set is not analyzed, because too few individuals are retained to fit the null model.\n\n", snpset);
			if(res==1) {
				fprintf(residf, "%s \n This SNP set is not analyzed, because too few individuals are retained to fit the null model.\n\n", snpset);
			}
			continue;
		}
		formatdata(n_fam, n_cov, n_snp, n_indiv, family, &y, &x, &G);
		nullResult = (struct NullResult *) malloc(sizeof(struct NullResult));
		SKAT_null (n_indiv, n_cov+1, res, y, x, lambda, U, nullResult);
		fprintf(paramf, "%s \nEstimate\tnullMLE\tSE_nullMLE\n", snpset);
		fprintf(paramf, "Heritability\t%lf\t%lf\n", (nullResult->error)[0][0], (nullResult->error)[0][1]);
		fprintf(paramf, "Additive_Var\t%lf\t%lf\n", (nullResult->error)[1][0], (nullResult->error)[1][1]);
		fprintf(paramf, "Error_Var\t%lf\t%lf\n", (nullResult->error)[2][0], (nullResult->error)[2][1]);
		fprintf(paramf, "Intercept\t%lf\t%lf\n", (nullResult->cov)[0][0], (nullResult->cov)[0][1]);
		for(i=1; i<n_cov+1; i++){
			fprintf(paramf, "Covariate_%d\t%lf\t%lf\n", i, (nullResult->cov)[i][0], (nullResult->cov)[i][1]);
		}
		fprintf(paramf, "\n");
		
		if(res==1){
			fprintf(residf, "%s \n", snpset);
			for(i=0; i<n_fam; i++){
				family_tmp = family+i;
				for(j=0; j<(family_tmp->n_indiv); j++){
					fprintf(residf, "%ld ", (family_tmp->indiv_included)[j]);
				}
			}
			fprintf(residf, "\n");
			for(i=0; i<n_indiv; i++){
				fprintf(residf, "%lf ", (nullResult->resid0)[i]);
			}
			fprintf(residf, "\n");
			for(i=0; i<n_indiv; i++){
				fprintf(residf, "%lf ", (nullResult->resid_std)[i]);
			}
			fprintf(residf, "\n\n");
		}

		MapRowXd U_map (U[0], n_indiv, n_indiv);
		MapRowXd G_map(G[0], n_indiv, n_snp);
		MapRowXd W_map(weight, n_snp, 1);
		MatrixXd UtGW_map = U_map.transpose() * G_map * W_map.asDiagonal(); // &UtGW_map(0) stores UtGW by columns			
		if(n_snp==1){
			printf("WARNING: Only 1 SNP present in SNP set %s! The MONSTER test is not performed.\n", snpset);
			fprintf(outf, "%s\t%d\t%d\tNA\tNA ", snpset, n_indiv, n_snp);
		} else{
			pvalue = SKAT (n_indiv, n_snp, n_grid, G, weight, nullResult, &UtGW_map(0), &rho_optim);
			fprintf(outf, "%s\t%d\t%d\t%.1lf\t%.12lf", snpset, n_indiv, n_snp, rho_optim, pvalue);
		}
		if(comparison==1){
			pvalue_burden = burdenTest(n_indiv, n_snp, G, weight, nullResult, &UtGW_map(0));
			pvalue_second = SecondOrderTest(n_indiv, n_snp, G, weight, nullResult, &UtGW_map(0));
			fprintf(outf, "\t%.12lf\t%.12lf", pvalue_burden, pvalue_second);
		}
		fprintf(outf, "\n");

		free(snpsetstart);
		free(weight);
		for(i=0; i<n_fam; i++){
			family_tmp = family+i;
			free(family_tmp->G);
			free_d2array(family_tmp->G_all, family_tmp->famsize, n_snp);
			free(family_tmp->indiv_included);
			free_d2array(family_tmp->kin_included, family_tmp->n_indiv, family_tmp->n_indiv);
		}

		if (erecyc==0) {
			free_d2array(U, n_indiv, n_indiv);
			free(lambda);
		}

		free(y);
		free_d2array(G, n_indiv, n_snp);
		free_d2array(x, n_indiv, n_cov+1);
		free(nullResult->resid);
		free_d2array(nullResult->Sigma_inv_half, n_indiv, n_indiv);
		free_d2array(nullResult->cov, n_cov+1, 2);
		if(res==1){
			free(nullResult->resid0);
			free(nullResult->resid_std);
		}
		free(nullResult);
		
	}
	fclose(snpf);
	fclose(outf);
	fclose(paramf);
	if(res==1) fclose(residf);
	
		
	for (i=0; i<n_fam; i++){
		free((family+i)->indiv);
		free((family+i)->pheno);
		free_d2array((family+i)->cov, (family+i)->famsize, n_cov);
		free_d2array((family+i)->kin, (family+i)->famsize, (family+i)->famsize);
		free((family+i)->missingrate);
	}
	free(family);
	if(erecyc == 1){
		free_d2array(U, n_indiv_all, n_indiv_all);
		free(lambda);
	}
	free(snpstart);

	for (i=0; i<n_allsnp; i++) free(allsnp[i]);
	free(allsnp);
	
	
	printf("\nProgram completed!\n");
	if(res==1 && efile ==1){
		printf("Output files: MONSTER.out\tMONSTER.param\tMONSTER.resid\n");
	} else if(res==1 && erecyc==1){
		printf("Output files: MONSTER.out\tMONSTER.param\tMONSTER.resid\tMONSTER.eig\n");
	} else if(res==1){
		printf("Output files: MONSTER.out\tMONSTER.param\tMONSTER.resid\n");
	} else if(efile==1){
		printf("Output files: MONSTER.out\tMONSTER.param\n");
	} else if(erecyc==1){
		printf("Output files: MONSTER.out\tMONSTER.param\tMONSTER.eig\n");
	} else{
		printf("Output files: MONSTER.out\tMONSTER.param\n");
	}

	return 0;	
}

int findSNP(int n_allsnp, char **allsnp, char snp[]){
	int index = 0;
	while (index < n_allsnp){
		if(strcmp(allsnp[index], snp)==0) return (index);
		index++;
	}
	return (-1);
}

void get_eigen(int n_fam, struct FAMILY *family, double **U, double *lambda){
	int i, famsize, j, index=0, k;
	struct FAMILY *family_tmp;
	for (i=0; i<n_fam; i++){
		family_tmp = family+i;
		famsize = family_tmp->famsize;
		MapRowXd phi_map ((family_tmp->kin)[0], famsize, famsize);
		SelfAdjointEigenSolver<MatrixXd> es;
		es.compute(phi_map);
		for (j=0; j<famsize; j++){
			lambda[index+j] = es.eigenvalues()(j);
			for(k=0; k<famsize; k++){
				U[index+j][index+k] = es.eigenvectors()(j, k);
			}
		}
		index += family_tmp->famsize;
	}
}


void get_eigen2(int n_fam, struct FAMILY *family, int *n_indiv, double ***U, double **lambda){
	int i, famsize, j, index=0, k, n_ind;
	struct FAMILY *family_tmp;
	*n_indiv = 0;
	for (i=0; i<n_fam; i++){
		*n_indiv += (family+i)->n_indiv;
	}
	(*U) = d2array(*n_indiv, *n_indiv);
	(*lambda) = d1array(*n_indiv);

	for (i=0; i<n_fam; i++){
		family_tmp = family+i;
		n_ind = family_tmp->n_indiv;
		if(n_ind==0)continue;
		MapRowXd phi_map ((family_tmp->kin_included)[0], n_ind, n_ind);
		SelfAdjointEigenSolver<MatrixXd> es;
		es.compute(phi_map);
		for (j=0; j<n_ind; j++){
			(*lambda)[index+j] = es.eigenvalues()(j);
			for(k=0; k<n_ind; k++){
				(*U)[index+j][index+k] = es.eigenvectors()(j, k);
			}
		}		
		index += n_ind;
	}	
}



void imputation(int n_fam, int n_snp, struct FAMILY *family, int BLUP, int BLUE, int AVE){
	void BLUP_impute(int n_fam, int snp_i, struct FAMILY *family);
	void BLUE_impute(int n_fam, int snp_i, struct FAMILY *family);
	void AVE_impute(int n_fam, int snp_i, struct FAMILY *family);
	
	int i;
	if(BLUP==1){
		for(i=0; i<n_snp; i++) BLUP_impute(n_fam, i, family);
	} else if(BLUE==1){
		for(i=0; i<n_snp; i++) BLUE_impute(n_fam, i, family);
	} else{
		for(i=0; i<n_snp; i++) AVE_impute(n_fam, i, family);
	}
}

void BLUP_impute(int n_fam, int snp_i, struct FAMILY *family){
	int i, n_indiv, n_indiv_nonmissing, j, flag, k, *nonmissing, flag2;
	struct FAMILY *family_tmp;
	double sum=0, sumG=0, freq;
	
	for (i=0; i<n_fam;i++){
		family_tmp = family+i;
		n_indiv = family_tmp->n_indiv;
		n_indiv_nonmissing=0;
		for(j=0; j<n_indiv; j++){
			if((family_tmp->G)[j][snp_i]!=-9) {
				n_indiv_nonmissing++;
				continue;
			}
		}			
		family_tmp->n_missing = family_tmp->n_indiv - n_indiv_nonmissing;
		if(n_indiv_nonmissing==0) {
			continue;
		}
		family_tmp->missing = i1array(family_tmp->n_missing);
		nonmissing = i1array(n_indiv_nonmissing);
		family_tmp->Gtmp = d1array(family_tmp->n_missing);
		flag = 0;
		flag2 = 0;
		for(j=0; j<n_indiv; j++){
			if((family_tmp->G)[j][snp_i]!=-9) {
				(nonmissing)[flag] = j;
				flag++;
			} else{
				(family_tmp->missing)[flag2] = j;
				flag2++;
			}
		}		
		MatrixXd phi_nonmissing(n_indiv_nonmissing, n_indiv_nonmissing);
		MatrixXd phi_MN(family_tmp->n_missing, n_indiv_nonmissing);
		VectorXd b = VectorXd::Constant(n_indiv_nonmissing, 1);
		VectorXd geno_nonmissing(n_indiv_nonmissing);
		for (j=0; j<n_indiv_nonmissing; j++){
			geno_nonmissing(j) = (family_tmp->G)[nonmissing[j]][snp_i];
			for(k=0; k<n_indiv_nonmissing; k++){
				phi_nonmissing(j, k) = (family_tmp->kin_included)[nonmissing[j]][nonmissing[k]];
			}
			for(k=0; k<(family_tmp->n_missing); k++){
				phi_MN(k, j) = (family_tmp->kin_included)[(family_tmp->missing)[k]][nonmissing[j]];
			}
		}
		LLT<MatrixXd> lltOfA(phi_nonmissing);
		lltOfA.solveInPlace(b);
		sum += b.sum();
		sumG += (b.transpose() * geno_nonmissing)(0);
		lltOfA.solveInPlace(geno_nonmissing);
		VectorXd tmp1 = phi_MN * geno_nonmissing;
		VectorXd tmp2 = phi_MN * b;
		for (j=0; j<(family_tmp->n_missing); j++) {
			(family_tmp->Gtmp)[j] = tmp2(j);
			(family_tmp->G)[(family_tmp->missing)[j]][snp_i] = tmp1(j);
		}
		free(nonmissing);
	}
	freq = sumG/sum;
	for(i=0; i<n_fam; i++){
		family_tmp = family+i;
		if(family_tmp->n_missing == family_tmp->n_indiv){
			for(j=0; j<family_tmp->n_indiv; j++){
				(family_tmp->G)[j][snp_i] = freq;
			}
		} else{
			for(j=0; j<family_tmp->n_missing; j++){
				(family_tmp->G)[(family_tmp->missing)[j]][snp_i] += freq * (1 - (family_tmp->Gtmp)[j]);
			}			
			free(family_tmp->Gtmp);
			free(family_tmp->missing);
		}
	}
		
}

void BLUE_impute(int n_fam, int snp_i, struct FAMILY *family){
	int i, n_indiv, n_indiv_nonmissing, j, *nonmissing, flag, k;
	struct FAMILY *family_tmp;
	double sum=0, sumG=0, freq;
	
	for (i=0; i<n_fam;i++){
		family_tmp = family+i;
		n_indiv = family_tmp->n_indiv;
		n_indiv_nonmissing=0;
		for(j=0; j<n_indiv; j++){
			if((family_tmp->G)[j][snp_i]!=-9) {
				n_indiv_nonmissing++;
				continue;
			}
		}
		nonmissing = i1array(n_indiv_nonmissing);
		flag = 0;
		for(j=0; j<n_indiv; j++){
			if((family_tmp->G)[j][snp_i]!=-9) {
				nonmissing[flag] = j;
				flag++;
				continue;
			}
		}		
		MatrixXd phi_nonmissing(n_indiv_nonmissing, n_indiv_nonmissing);
		VectorXd b = VectorXd::Constant(n_indiv_nonmissing, 1);
		MatrixXd geno_nonmissing(n_indiv_nonmissing, 1);
		for (j=0; j<n_indiv_nonmissing; j++){
			geno_nonmissing(j,0) = (family_tmp->G)[nonmissing[j]][snp_i];
			for(k=0; k<n_indiv_nonmissing; k++){
				phi_nonmissing(j, k) = (family_tmp->kin_included)[nonmissing[j]][nonmissing[k]];
			}
		}
		LLT<MatrixXd> lltOfA(phi_nonmissing);
		lltOfA.solveInPlace(b);
		sum += b.sum();
		sumG += (b.transpose() * geno_nonmissing)(0,0);
		free(nonmissing);
	}
	freq = sumG/sum;
	for (i=0; i<n_fam; i++){
		family_tmp = family+i;
		for(j=0; j<family_tmp->n_indiv; j++){
			if((family_tmp->G)[j][snp_i]==-9) (family_tmp->G)[j][snp_i]=freq;
		}
	}	
	
}

void AVE_impute(int n_fam, int snp_i, struct FAMILY *family){
	double geno=0, ave;
	int i, j, n_nonmissing=0;
	struct FAMILY *family_tmp;
	for(i=0; i<n_fam; i++){
		family_tmp = family+i;
		for(j=0; j<family_tmp->n_indiv; j++){
			if((family_tmp->G)[j][snp_i] != -9){
				n_nonmissing++;
				geno += (family_tmp->G)[j][snp_i];
			}
		}
	}
	ave = geno/n_nonmissing;
	for(i=0; i<n_fam; i++){
		family_tmp = family+i;
		for(j=0; j<family_tmp->n_indiv; j++){
			if((family_tmp->G)[j][snp_i] == -9){
				(family_tmp->G)[j][snp_i] = ave;
			}
		}
	}	

}


void formatdata(int n_fam, int n_cov, int n_snp, int n_indiv, struct FAMILY *family, double **y, double ***x, double ***G){
	int i, j, flag=0, k;
	struct FAMILY *family_tmp;
	(*y) = d1array(n_indiv);
	(*x) = d2array(n_indiv, n_cov+1);
	(*G) = d2array(n_indiv, n_snp);
	for(i=0; i<n_fam; i++){
		family_tmp = family+i;
		for(j=0; j<family_tmp->n_indiv; j++){
			(*y)[flag] = (family_tmp->pheno)[j];
			(*x)[flag][0] = 1;
			for (k=0; k<n_snp; k++){
				(*G)[flag][k] = (family_tmp->G)[j][k];
			}
			for(k=0; k<n_cov; k++){
				(*x)[flag][k+1] = (family_tmp->cov)[j][k];
			}
			flag++;
		}
	}
}


void formatdata_null(int n_fam, int n_cov, int n_indiv, struct FAMILY *family, double **y, double ***x){
	int i, j, flag=0, k;
	struct FAMILY *family_tmp;
	(*y) = d1array(n_indiv);
	(*x) = d2array(n_indiv, n_cov+1);
	for(i=0; i<n_fam; i++){
		family_tmp = family+i;
		for(j=0; j<family_tmp->famsize; j++){
			(*y)[flag] = (family_tmp->pheno)[j];
			(*x)[flag][0] = 1;
			for(k=0; k<n_cov; k++){
				(*x)[flag][k+1] = (family_tmp->cov)[j][k];
			}
			flag++;
		}
	}
}





