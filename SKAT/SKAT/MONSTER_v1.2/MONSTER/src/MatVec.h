#ifndef _MatVec_H_
#define _MatVec_H_

# include <stdlib.h>
#include <stdio.h>
# include <math.h>


/************************************************************************
******** version 5/20/2012, new functions: ******************************
****************** matrixIndvLU, lubksb, ludcmp *************************
****************** matrixInvCholesky, cholsl, choldc ********************
**************************************************************************/


# define FREE_ARG char*


double *d1array(int n)
// allocate a double vector with dimention n
{
	double *x;
	int i=0;
	x=(double *)malloc(n*sizeof(double));
	for (i=0;i<n;i++){
		x[i]=0;
	}
	return x;
}

int *i1array(int n)
// allocate an int vector with dimention n
{
	int *x;
	int i;
	x=(int *)malloc(n*sizeof(int));
	for (i=0;i<n;i++){
		x[i]=0;
	}
	return x;
}

void showD1array(double *x, int n)
{
	for (int i=0; i<n;i++) printf("%lf\t",x[i]);
	printf("\n");
}

void showI1array(int *x, int n)
{
	for (int i=0; i<n;i++) printf("%d\t",x[i]);
	printf("\n");
}

double **d2array(int m,int n)
// allocate a double matrix with dimensions m and n
{
	int i, j;
	double **x;
	x=(double **) malloc(m*sizeof(double*));
	x[0]=(double *) malloc(m*n*sizeof(double));
	for(i=1;i<m;i++) x[i]=x[i-1]+n;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			x[i][j]=0;
		}
}
return x;
}

int **i2array(int m,int n)
// allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]
{
	int i, j;
	int **x;
	x=(int **) malloc(m*sizeof(int*));
	x[0]=(int *) malloc(m*n*sizeof(double));
	for(i=1;i<m;i++) x[i]=x[i-1]+n;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			x[i][j]=0;
		}
}
return x;
}

void showD2array(double **x, int m, int n)
{
	for (int i=0; i<m;i++) {
		for (int j=0; j<n; j++){
			printf("%lf\t",x[i][j]);
		}
		printf("\n");
	}
}

void showI2array(int **x, int m, int n)
{
	for (int i=0; i<m;i++) {
		for (int j=0; j<n; j++){
			printf("%d\t",x[i][j]);
		}
		printf("\n");
	}
}

void free_d2array(double **x, int m, int n)
// free a double matrix with dimensions m and n
{
free((FREE_ARG) (x[0]));
free((FREE_ARG) (x));
}

void free_i2array(int **x, int m, int n)
// free a double matrix with dimensions m and n
{
free((FREE_ARG) (x[0]));
free((FREE_ARG) (x));
}



double **constD2arrayGen (double c, int m, int n)
{
	double **result;
	int i, j;
	result=d2array(m,n);
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++){
			result [i] [j] = c;
		}
	}
	return result;
}


double *tranToD1array (double **mat, int m, bool type)
//transform matrix to vector
// type=0 if mat is a row
// type=1 if mat is a column
{
	int i;
	double *result=d1array(m);
	if (!type) {
		for (i=0; i<m; i++) result[i]=mat[0][i];
	}
	else {
		for (i=0; i<m; i++) result[i]=mat[i][0];
	}
	return result;
}

double d1arrayMax(double *A, int n)
{
	double max=A[0];
	int i;
	for (i=0;i<n;i++){
		if (A[i]>max) max=A[i];
	}
	return max;
}

double d1arrayMin(double *A, int n)
{
	double min=A[0];
	int i;
	for (i=0;i<n;i++){
		if (A[i]<min) min=A[i];
	}
	return min;
}

int d1arrayMinIndex(double *A, int n)
{
	double min=A[0];
	int i,flag=1;
	for (i=0;i< n;i++){
		if (A[i]<min) {
			flag=i;
			min=A[i];
		}
	}
	return flag;
}


double d2arrayMax(double **A, int m, int n)
{
	double max=A[0][0];
	int i,j;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			if (A[i][j]>max) max=A[i][j];
		}
	}
	return max;
}

int *d2arrayMaxIndex(double **A, int m, int n, int *flag)
{
	double max=A[0][0];
	int i,j;
	for (i=0;i< m;i++){
		for (j=0; j<n;j++){
			if (A[i][j]>max) {
				flag[0]=i;
				flag[1]=j;
				max=A[i][j];
			}
		}
	}
	return flag;
}


double d1arraySum (double *A, int n)
{
	double sum=0;
	int i;
	for (i=0;i<n;i++) sum+=A[i];
	return sum;
}

int i1arraySum (int *A, int n)
{
	int sum=0;
	int i;
	for (i=0;i<n;i++) sum+=A[i];
	return sum;
}


double d2arraySum(double **A, int m, int n)
{
	double sum=0;
	int i,j;
	for (i=0;i<m;i++){
		for (j=0; j<n; j++){
			sum+=A[i][j];
		}
	}
	return sum;
}

double d1arrayMean (double *A, int n)
{
	double sum=0;
	int i;
	for (i=0;i<n;i++) sum+=A[i];
	return sum/n;
}

double d2arrayMean(double **A, int m, int n)
{
	double sum=0;
	int i,j;
	for (i=0;i<m;i++){
		for (j=0; j<n; j++){
			sum+=A[i][j];
		}
	}
	return sum/(m*n);
}

double* d2arrayRowSums(double **A, int m, int n)
{
	double *Sum;
	int i,j;
	Sum=d1array(m);
	for (i=0; i<m; i++){
		Sum[i]=0;
		for (j=0; j<n; j++){
			Sum[i]+=A[i][j];
		}
	}
	return Sum;
}

double* d2arrayColSums(double **A, int m, int n)
{
	double *Sum;
	int i,j;
	Sum=d1array(n);
	for (j=0; j<n; j++){
		Sum[j]=0;
		for (i=0; i<m; i++){
			Sum[j]+=A[i][j];
		}
	}
	return Sum;
}

double* d2arrayRowMeans(double **A, int m, int n)
{
	double *Sum;
	int i,j;
	Sum=d1array(m);
	for (i=0; i<m; i++){
		Sum[i]=0;
		for (j=0; j<n; j++){
			Sum[i]+=A[i][j];
		}
		Sum[i]=Sum[i]/n;
	}
	return Sum;
}

double* d2arrayColMeans(double **A, int m, int n)
{
	double *Sum;
	int i,j;
	Sum=d1array(n);
	for (j=0; j<n; j++){
		Sum[j]=0;
		for (i=0; i<m; i++){
			Sum[j]+=A[i][j];
		}
		Sum[j]=Sum[j]/m;
	}
	return Sum;
}

double ** ScalarProd(double s, void *A, int m, int n, bool typeA)
{
	int i,j;
	double **sA=d2array(m,n);
	if (typeA){
		for (i=0;i<m;i++){
			for (j=0;j<n;j++){
				sA[i][j]=s*((double**)A)[i][j];
			}
		}
	} else {
		if (m==1){
			for (j=0;j<n;j++)sA[0][j]=s*((double*)A)[j];
		} else {
			for (i=0; i<m;i++)sA[i][0]=s*((double*)A)[i];
		}
	}
	return sA;
}

double **darrayDiff(void *A, void *B, int dimI, int dimJ, bool typeA, bool typeB)
// tapeA=1 if A is a double ** matrix
//              0 if A is a double * vector
// same for typeB
{
	double **AB;
	int i, j;
	AB=d2array(dimI,dimJ);

	if ( typeA && typeB ) {
		for (i=0; i<dimI; i++){
			for (j=0; j<dimJ; j++){
				AB [i][j] =  ((double **)A) [i][j] - ((double **)B) [i][j];
			}
		}
		return AB;
	}
	if ( typeA && !( typeB) ) {
		if (dimJ==1){
			for (i=0; i<dimI; i++){
				AB[i][0] = ((double **)A) [i][0] - ((double *) B) [i];
			}
		} else {
			for (j=0; j<dimJ; j++){
				AB[0][j] += ((double **)A) [0][j] - ((double *) B) [j];
			}
		}
		return AB;
	}

	if ( !(typeA) && typeB ) {
		if (dimJ==1){
			for (i=0; i<dimI; i++){
				AB[i][0] += ((double *) A) [i] - ((double **)B) [i][0];
			}
		} else {
			for (j=0; j<dimJ; j++){
				AB[0][j] += ((double *) A) [j] - ((double **)B) [0][j];
			}
		}
		return AB;
	}

	if ( !(typeA) && !(typeB) ){
		if (dimJ==1){
			for (i=0; i<dimI; i++){
				AB[i][0] += ((double *) A) [i] - ((double *) B) [i];
			}
		} else {
			for (j=0; j<dimJ; j++) AB[0][j] += ((double *) A) [j] - ((double *) B) [j];
						}
		return AB;
	}




}



double **  darrayMult (void *A, void *B, int dimI, int dimJ, int dimK, bool typeA, bool typeB)
// tapeA=1 if A is a double ** matrix
//              0 if A is a double * vector
// same for typeB
{
	double **AB;
	int i, j, k;
	AB=d2array(dimI,dimK);

	if ( typeA && typeB ) {
		for (i=0; i<dimI; i++){
			for (j=0; j<dimJ; j++){
				for (k=0; k<dimK; k++){
					AB [i][k] +=  ((double **)A) [i][j] * ((double **)B) [j][k];
				}
			}
		}
		return AB;
	}
	if ( typeA && !( typeB) ) {
		if (dimJ==1){
			for (i=0; i<dimI; i++){
				for (k=0; k<dimK; k++){
					AB[i][k] += ((double **)A) [i][0] * ((double *) B) [k];
				}
			}
		} else {
			for (i=0; i<dimI; i++){
				for (j=0; j<dimJ; j++){
					AB[i][0] += ((double **)A) [i][j] * ((double *) B) [j];
				}
			}
		}
		return AB;
	}

	if ( !(typeA) && typeB ) {
		if (dimJ==1){
			for (i=0; i<dimI; i++){
				for (k=0; k<dimK; k++){
					AB[i][k] += ((double *) A) [i] * ((double **)B) [0][k];
				}
			}
		} else {
			for (j=0; j<dimJ; j++){
				for (k=0; k<dimK; k++){
					AB[0][k] += ((double *) A) [j] * ((double **)B) [j][k];
				}
			}
		}
		return AB;
	}

	if ( !(typeA) && !(typeB) ){
		if (dimJ==1){
			for (i=0; i<dimI; i++){
				for (k=0; k<dimK; k++){
					AB[i][k] += ((double *) A) [i] * ((double *) B) [k];
				}
			}
		} else {
			for (j=0; j<dimJ; j++) AB[0][0] += ((double *) A) [j] * ((double *) B) [j];
						}
		return AB;
	}




}

double ** arrayTranspose(double **x, int m, int n)
{
	double **result=d2array(n,m);
	int i,j;
	for (i=0;i<n;i++){
		for (j=0;j<m;j++){
			result[i][j]=x[j][i];
		}
	}
	return result;
}

#define TINY 1.0e-20
void matrixInvLU(double **a, int n, double **b){
//find matrix inverse through LU decomposition
// a is matrix to be inversed, n is its dimension
// result stored in b
    void lubksb(double **a, int n, int *indx, double b[]);
    void ludcmp(double **a, int n, int *indx, double *d);
    int i,j, *indx;
    double **aa,d0=0,*d=&d0;
    indx=i1array(n);
    aa=d2array(n,n);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            aa[i][j]=a[j][i];
            }
    }

    ludcmp(aa,n,indx,d);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            b[i][j]=0;
        }
        b[i][i]=1;
        lubksb(aa,n,indx,b[i]);
	}
    free_d2array(aa,n,n);
    free(indx);
}

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii-1;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
	vv=d1array(n);
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) printf("ERROR: Singular matrix in routine ludcmp!\n");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != (n-1)) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	free(vv);
}


void matrixInvCholesky(double **a, int n, double **ans){
//find the inverse of a positive symmetric definite through Cholesky decomposition
// a is matrix to be inversed, n is its dimension
// result stored in b
    void cholsl(double **a, int n, double p[], double b[], double x[]);
    void choldc(double **a, int n, double p[]);
    int i,j, *indx;
    double **aa, *p,*b;
    p=d1array(n);
    aa=d2array(n,n);
    b=d1array(n);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            aa[i][j]=a[i][j];
            }
    }
    choldc(aa,n,p);
    for (i=0;i<n;i++){
        for (j=0;j<n;j++) b[j]=0;
        b[i]=1;
        cholsl(aa,n,p,b,ans[i]);
	}
    free_d2array(aa,n,n);
    free(p);
    free(b);
}
void cholsl(double **a, int n, double p[], double b[], double x[])
/*Solves the set of n linear equations A ， x = b, where a is a positive-definite symmetric matrix.
a[0..n-1][0..n-1] and p[0..n-1] are input as the output of the routine choldc. Only the lower
subdiagonal portion of a is accessed. b[0..n-1] is input as the right-hand side vector. The
solution vector is returned in x[0..n-1]. a, n, and p are not modified and can be left in place
for successive calls with different right-hand sides b. b is not modified unless you identify b and
x in the calling sequence, which is allowed.
*/
{
    void choldc(double **a, int n, double p[]);
    int i,k;
    double sum;
    for (i=0;i<n;i++) { //Solve L ， y = b, storing y in x.
        for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
        x[i]=sum/p[i];
    }
    for (i=n-1;i>=0;i--) { //Solve LT ， x = y.
        for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
        x[i]=sum/p[i];
    }
}


void choldc(double **a, int n, double p[])
/*Given a positive-definite symmetric matrix a[0..n-1][0..n-1], this routine constructs its Cholesky
decomposition, A = L ， LT . On input, only the upper triangle of a need be given; it is not
modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
elements which are returned in p[1..n].
*/
{
    int i,j,k;
    double sum;
    for (i=0;i<n;i++) {
        for (j=i;j<n;j++) {
            for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
            if (i == j) {
                if (sum <= 0.0) //a, with rounding errors, is not positive definite.
                    printf("Choldc failed\n");
                p[i]=sqrt(sum);
            } else a[j][i]=sum/p[i];
        }
    }
}



#endif

