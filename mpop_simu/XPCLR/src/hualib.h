
#ifndef __hclib__
#define __hclib__


#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define NR_END 1
#define FREE_ARG char*
#define IA  16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define CLK_TCK CLOCKS_PER_SEC
#define DBL_MAX 1.79769313486231470e+308
#define PI 3.141592654

#define ERROR(message) fprintf(stderr,message),fprintf(stderr,"\n"), exit(1)
#define WARNNING(message) fprintf(stderr,message),fprintf(stderr,"\n")

#define min(x,y) (((x) <(y)) ? x: y)
#define max(x,y) (((x) >(y)) ? x: y)
#define MIN(x,y) ((x)<(y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y) ? (x) : (y) )
//ONLY for comparing char and int variables, and returns "1" if equal. 
#define EQUAL(x,y) ((x)==(y) ? 1 : 0 ) 
#define GEd(x,y) ((x)>=(y) ? 1.0 : 0.0)
#define GEi(x,y) ((x)>=(y) ? 1 : 0)

 void nrerror(char error_text[]);

 
 /******************************************************************************
   subroutines for dynamic memory allocation.
 *******************************************************************************/
double *vector(long nl, long nh);
 double *dvector(long nl, long nh);
 int *ivector(long nl, long nh);
 char *cvector(long nl, long nh);
 double **matrix(long nrl, long nrh, long ncl, long nch);
 void matrixvector(double *v1, double **m, double *v2, int n);
 void zerovec(double *v, int n);
 double **dmatrix(long nrl, long nrh, long ncl, long nch);
 char **cmatrix(long nrl, long nrh, long ncl, long nch);
 int **imatrix(long nrl, long nrh, long ncl, long nch);
 int ***i3matrix(long nrl, long nrh, long ncl, long nch, long n3l, long n3h);
 double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
 double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);

 double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
 void free_vector(double *v, long nl, long nh);
 void free_ivector(int *v, long nl, long nh);
 void free_cvector(char *v, long nl, long nh);
// void free_cvector(unsigned char *v, long nl, long nh);
 void free_lvector(unsigned long *v, long nl, long nh);
 void free_dvector(double *v, long nl, long nh);
 void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
 void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
 void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
 void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
 void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
 void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
 
 /****************************************************************************
   subroutines for etc.
 ******************************************************************************/
   
 double log_pull(double logp, double logq);
 char *utoa(unsigned value, char *digits, int base);
 char *itoa(int value, char *digits, int base);


 /*****************************************************************************
         subroutines for sorting.
 *******************************************************************************/  
 void qs(long *items,long *rank, int left, int right);
int mnmial(int n,int nclass,double p[],int rv[]);
int ordran(int n,double pbuf[]);
int order(int n,double pbuf[]);
int ranvec(int n,double pbuf[]);


/*******************************************************************************
        dynamic allocation memories.
*******************************************************************************/

template < typename T >
T **Allocate2DArray( int nRows, int nCols);
template < typename T >
void Free2DArray(T** Array);


#ifndef __hua_stat__
#define __hua_stat__


/*****************************************************************************
   subroutines for generating random variables.
******************************************************************************/   

 double ran1(int *idum);
 double ran1();
double binopdf(int n,int k,double p);
double factln(int n);
double bico(int n,int k);

 double exprnd(double rate, int *idum);
  /*generate binomial rv. */
 double bnldev(double pp, int n, int *idum);
 double betadev(double a, double b, int *idum);
 double gamdev(int ia, int *idum);
 double gamma(double shape, double scale, int *idum);
/*Normal(0,1)*/
 double gasdev(int *idum);
double gasdev(double m,double v);

int poisso(double u);
 double gammln(double xx);

void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammp(double a, double x);
double erff(double x);
double gaussiancdf(double x, double m, double s);


double qtrap(double (*func)(double),double a,double b);
double trapzd(double (*func)(double),double a,double b,int n);
double qromb(double (*func)(double),double a,double b);
void polint(double xa[],double ya[],int n,double x,double *y,double *dy);




/***************************************************************************
                  fundermental statistical subroutines.
****************************************************************************/
 void moment(double data[], int n, double *ave, double *adev, double *sdev,
              double *var, double *skew, double *curt);


/****************************************************************************
                  statistical tests.
*****************************************************************************/
/*     chi-square test   */
 void chsone(float bins[], float ebins[], int nbins, int knstrn, float *df,
	float *chsq, float *prob);


/****                Kolmogorov-Smirnov Test               ****/
void ksone(float data[], unsigned long n, float (*func)(float), float *d,
	float *prob);


#endif


#endif
