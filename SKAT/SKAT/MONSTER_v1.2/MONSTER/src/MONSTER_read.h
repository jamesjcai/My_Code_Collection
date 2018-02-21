#ifndef _MONSTER_READ_H_
#define _MONSTER_READ_H_

#define MAXCOV 100
#define MAXSNPLEN 30
#define MAXINDIV 10000

struct FAMILY {
	long int famID;
	int famsize, n_indiv, *missing, n_missing;
	long int *indiv, *indiv_included;
	double **kin, *missingrate, **kin_included, **eigvector, *eigvalue, *Gtmp;
	double *pheno, **cov, **G_all, **G; 
};

void count_ped (char *pedfile, int *n_indiv_all, int *n_cov, int *n_fam){
	FILE *ped;
	int length=20*MAXCOV+150, n=0, tmp, i;	
	long int famID_prev, famID, ID_tmp;
	double cov_tmp, pheno_tmp;
	char line[20*MAXCOV+150];	
	
	*n_indiv_all = 0;
	*n_cov = 0;
	*n_fam = 0;
	
	if ((ped = fopen(pedfile, "r")) == NULL ){
		printf("ERROR: Failure to open pedigree and phenotype data file!\n");
		exit(1);
	}
    fgets(line, length, ped);
	sscanf(line,"%ld %ld %ld %ld %d %[^\n]", &famID_prev, &ID_tmp, &ID_tmp, &ID_tmp, &tmp, line);
    n = sscanf(line,"%lf %[^\n]", &pheno_tmp, line);
    if (n<1) {	
		printf("WARNING: No covariate in pedigree and phenotype data file!\n");
    } else {
		while (n==2){
			n=sscanf(line,"%lf %[^\n]", &cov_tmp, line);
			if (n>=1)(*n_cov)++;		
		}
    }
	(*n_fam)++;
	(*n_indiv_all)++;
	while (!feof(ped)){
		if(fscanf(ped, "%ld %ld %ld %ld %d %lf", &famID, &ID_tmp, &ID_tmp, &ID_tmp, &tmp, &pheno_tmp) != 6){
			break;
		}
		(*n_indiv_all)++;
		if(famID != famID_prev) {
			(*n_fam)++;
			famID_prev = famID;
		}
		fgets(line, length, ped);
	}	
	fclose(ped);	
}

int read_ped(char *pedfile, int n_indiv_all, int n_fam, int n_cov, struct FAMILY *family){
	FILE *ped;	
	int i, tmp, fam_index=-1, j, k, maxfamsize=0;
	long int famID, famID_prev=-1, indivID, ID_tmp;
	double pheno_tmp, cov_tmp;
	struct FAMILY *family_tmp;
	ped = fopen(pedfile, "r");
	for (i=0; i<n_fam; i++) (family+i)->famsize = 0;	
	while (!feof(ped)){
		if(fscanf(ped, "%ld %ld %ld %ld %d %lf", &famID, &indivID, &ID_tmp, &ID_tmp, &tmp, &pheno_tmp)!=6) break;
		if(famID != famID_prev){
			fam_index++;
			(family+fam_index) -> famsize = 1;
			(family+fam_index) -> famID = famID;
			famID_prev = famID;
		} else{
			(family+fam_index) -> famsize ++;
		}
		for (i=0; i<n_cov; i++){
			fscanf(ped, "%lf", &cov_tmp);
		}
	}
	fclose(ped);
	
	ped = fopen(pedfile, "r");
	for (i=0; i<n_fam; i++){
		family_tmp = family + i;
		family_tmp -> pheno = d1array(family_tmp->famsize);
		family_tmp -> cov = d2array(family_tmp->famsize, n_cov+1);		
		family_tmp->indiv = (long int*)malloc(family_tmp->famsize * sizeof(long int));
		family_tmp->missingrate = d1array(family_tmp->famsize);
		for (j=0; j<(family_tmp->famsize); j++){
			fscanf(ped, "%ld %ld %ld %ld %d %lf", 
				&famID, &((family_tmp->indiv)[j]), &ID_tmp, &ID_tmp, &tmp, &((family_tmp->pheno)[j]));
			for (k=0; k<n_cov; k++){
				fscanf(ped, "%lf", &((family_tmp->cov)[j][k]));
			}
		}
		if(family_tmp->famsize > maxfamsize) maxfamsize = family_tmp->famsize;
	}	
	fclose(ped);
	
	return maxfamsize;
	
}

int findInd(int n_indiv, long int *indivlist, long int indiv);
void read_kin(char *kinfile, int n_fam, struct FAMILY *family){

	int fam_i=-1, fam_old=-1, i, index1, index2, index1old=-1, index2old=-1;
	long int famID, famID_prev=-1, indiv1, indiv2, indiv1old=-1, indiv2old=-1;
	double coef_tmp;
	struct FAMILY *family_tmp;
	FILE *kin;
	if( (kin = fopen(kinfile, "r"))==NULL) {
		printf("ERROR: Failure to open the kinship coefficient file!\n");
		exit(1);
	}
	for (i=0; i<n_fam; i++){
		family_tmp = family+i;
		family_tmp->kin = d2array(family_tmp->famsize, family_tmp->famsize);
	}
	while (!feof(kin)){
		if(fscanf(kin, "%ld %ld %ld %lf", &famID, &indiv1, &indiv2, &coef_tmp)!=4) break;
		if(famID == famID_prev){
			if(indiv1 == indiv1old){
				index1 = index1old;
			} else if(indiv1 == indiv2old){
				index1 = index2old;
			} else {
				index1 = findInd(family_tmp->famsize, family_tmp->indiv, indiv1);
			}
			if(indiv2 ==indiv2old){
				index2 = index2old;
			} else if(indiv2 == indiv1old){
				index2 = index1old;
			} else {
				index2 = findInd(family_tmp->famsize, family_tmp->indiv, indiv2);
			}
			
		} else{
			fam_i++;
			family_tmp = family+fam_i;
			index1 = findInd(family_tmp->famsize, family_tmp->indiv, indiv1);
			index2 = findInd(family_tmp->famsize, family_tmp->indiv, indiv2);
		}
		if(index1 != index2){
			(family_tmp->kin)[index1][index2] = coef_tmp * 2;
			(family_tmp->kin)[index2][index1] = coef_tmp * 2;
		} else{
			if(coef_tmp>0.1){
				printf("WARNING: Inbreeding coefficient for individual %ld is %lf > 0.1. This is rare for most human samples. ", indiv1, coef_tmp);
				printf("User should make sure the kinship input file is correctly prepared. In particular, the file should contain inbreeding coefficients, not self-kinship coefficients.\n\n");
			}
			(family_tmp->kin)[index1][index2] = coef_tmp + 1;
		}
		famID_prev= famID;
		index1old = index1;
		index2old = index2;
		indiv1old = indiv1;
		indiv2old = indiv2;
	}
}

int findInd(int n_indiv, long int *indivlist, long int indiv){
	int index = 0;
	while (index < n_indiv){
		if(indivlist[index] == indiv) return (index);
		index++;
	}
	printf("ERROR: Individual %ld is present in the kinship file, but not in the phenotype and pedigree file!\n", indiv);
	exit(1);
}


int read_eig(char *eigfile, int n_indiv_all, double **U, double *lambda)	{
	FILE *eig;
	int i, j;
	if ((eig = fopen(eigfile, "r")) == NULL ){
		return 0;
	}			
	for (i=0; i<n_indiv_all; i++){
		fread(&(lambda[i]), sizeof(double), 1, eig);
	}
	for (i=0; i<n_indiv_all; i++){
		for (j=0; j<n_indiv_all; j++){
			fread(&(U[i][j]), sizeof(double), 1, eig);		
		}
	}
	fclose(eig);
	return 1;
}


int count_allsnp(char *genofile){
	FILE *geno;
	int n_allsnp=-1, length=MAXINDIV*10+MAXSNPLEN;
	char line[MAXINDIV*10+MAXSNPLEN];
	
	if( (geno = fopen(genofile, "r"))==NULL) {
		printf("ERROR: Failure to open the genotype data file!\n");
		exit(1);
	}	
	while (!feof(geno)){
	    if (fgets(line, length, geno)==NULL) break;
		n_allsnp++;
	}
	fclose(geno);
	return n_allsnp;
}

void read_allsnp(char *genofile, int n_fam, struct FAMILY *family, char **allsnp, int *snpstart){
	FILE *geno;
	int length=MAXINDIV*10+MAXSNPLEN, i, len, j, k;
	char line[MAXINDIV*10+MAXSNPLEN];
	long int indivID;
	struct FAMILY *family_tmp;
	
	geno = fopen(genofile, "r");
	fgets(line, length, geno);	
	len = strlen(line);
	snpstart[i] = len;
	
	sscanf(line, "%ld %[^\n]", &indivID, line);
	for(j=0; j<n_fam; j++){
		family_tmp = family+j;
		for(k=0; k<(family_tmp->famsize); k++){
			sscanf(line, "%ld %[^\n]", &indivID, line);
			if(indivID != (family_tmp->indiv)[k]) {
				printf("ERROR: Ordering of individuals in the genotype data file does not match that in the pedigree file!\n");
				exit(1);
			}
		}	
	}	

	i=0;
	while (!feof(geno)){
	    if (fgets(line, length, geno)==NULL) break;
		snpstart[i] = len;
		len += strlen(line);
		sscanf(line,"%s %[^\n]", allsnp[i], line);
		i++;
	}
	fclose(geno);
}

void read_geno(char *genofile, int n_fam, int n_snp, int n_indiv_all, double mcut, int *snpsetstart, struct FAMILY *family){
	FILE *geno;
	int i, j, k, indiv, indiv2, l;
	long int indivID;
	char SNP[MAXSNPLEN];
	struct FAMILY *family_tmp;
	double geno_tmp, snp_missing;

	geno = fopen(genofile, "r");
	for(i=0; i<n_snp; i++){
		snp_missing = 0;
		fseek(geno, snpsetstart[i], SEEK_SET);
		fscanf(geno, "%s", SNP);
		for (j=0; j<n_fam; j++){
			family_tmp = family+j;
			for(k=0; k<(family_tmp->famsize); k++){	
				fscanf(geno, "%lf", &geno_tmp);
				(family_tmp->G_all)[k][i] = geno_tmp;
				if(geno_tmp == -9) {
					(family_tmp->missingrate)[k] += 1;
					snp_missing++;
				}
			}
		}
		if(snp_missing>.5*n_indiv_all){
			printf("WARNING: For SNP %s, over 50%% of the individuals have missing genotypes!\n", SNP);
		}
	}
	fclose(geno);	
	for(j=0; j<n_fam; j++){
		family_tmp = family+j;
		family_tmp->n_indiv = 0;
		for(k=0; k<(family_tmp->famsize); k++){	
			(family_tmp->missingrate)[k] /= n_snp;
			if((family_tmp->missingrate)[k] <= mcut) family_tmp->n_indiv ++;
		}
		family_tmp->G = (double **)malloc(family_tmp->n_indiv * sizeof(double*));//**
		family_tmp->indiv_included = (long int*)malloc(family_tmp->n_indiv * sizeof(long int));
		family_tmp->kin_included = d2array(family_tmp->n_indiv, family_tmp->n_indiv);
		indiv = 0;
		for(k=0; k<(family_tmp->famsize); k++){	
			if((family_tmp->missingrate)[k] <= mcut){
				(family_tmp->G)[indiv] = (family_tmp->G_all)[k];//**
				(family_tmp->indiv_included)[indiv] = (family_tmp->indiv)[k];
				indiv2 = 0;
				for(l=0; l<(family_tmp->famsize); l++){	
					if((family_tmp->missingrate)[l] <= mcut){
						(family_tmp->kin_included)[indiv][indiv2] = (family_tmp->kin)[k][l];
						indiv2++;
					}
				}
				indiv++;
			}
		}
	}

}


#endif