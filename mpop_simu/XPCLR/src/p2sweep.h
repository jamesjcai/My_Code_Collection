
#ifndef __p2sweeplib__
#define __p2sweeplib__

#include "hualib.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#ifndef __MINMAX__
#define __MINMAX__

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

#endif //__MINMAX__


#define MINSITES 500

void getCount(int n1, int n2,int rawData[][10000],double *f1,double *f2,int *x1,int *x2,double *position, double *rawPosition,int nsegsites, int *nnsegsites,double rho);
void getCount2(int n1, int n2,double* rawq1, double* rawq2, double* f1, double *f2, int *x1, int *x2,double *rawposition, double *position,int nsegsites, int *nnsegsites);
double getW(double* q2, double* q1, int nsegsites);
double estimatec_s(double p0, double c, double s);
double func10(double t);
double func1(double t);
double func2(double t);
double func3(double t);
double logPNeutrality(double w,double p0,double *q2, int *x1,int n1,int nnsegsites);
double logPNeutrality2(double w,double p0,double *q2, int *x1,int n1,int *snpList,int SnpListLength);
double prob_sweep_fun(int n1, int x1,double q2,double c,double s, double p0,double w);
double logPSelection(double w,double p0,int *x1, double *q2, double *pos,double *sValue, int sLength,int n1, int indicator,int nnsegsites,double*logPLocation,double*maxS);
double logPSelection2(double w,double p0,int *x1, double *q2, double *pos,double *sValue, int sLength,int n1, int indicator,int nnsegsites,double* logPLocation, double* maxS, int *snpList,int SnpListLength,double gridsize,int gridnumber,double gridstart,double totalGD);
double logPNeutralityScan(double w,double p0,double *q2,int *x1,int *n1,int startP, int endP);
double logPSelectionScan(double w,double p0,int *x1,double *q2,double *pos,int *n1,int mutPos,int startP,int endP,double *maxS,double *sValue,int sLength);
double logPNeutralityScan2(double w,double p0,double *q2, int *x1,int *n1,int *snpList, int SnpListLength);
double logPSelectionScan2(double w,double p0,int *x1, double *q2, double *pos,int *n1, double mutGpos,double *maxS,double *sValue,int sLength, int *snpList, int SnpListLength);
int getSnpWindow(double *pos,int snpN,int *snpList, int *snpListLength, int windowSnpNumber,double mutGpos, double geneticWindowSize);
int getSubsetSnpWindow(int windowSnpNumber, double geneticWindowSize,int *snpList, int *snpListLength, int *idum);
int getSubsetSnpWindow_ms(int windowSnpNumber,int *snpList,int *snpListLength,int *idum);
int getSnpWindow_ms(int snpN, int *snpList, int *snpListLength, int windowSnpNumber);
double fst(double *q2,double *q1,int n2, int n1,int nnsegsites);
double logPNeutralityScan3(double w,double p0,double *q2, int *x1,int *n1,int *snpList, int SnpListLength,double *wt);
double logPSelectionScan3(double w,double p0,int *x1, double *q2, double *pos,int *n1, double mutGpos,double *maxS,double *sValue,int sLength, int *snpList, int SnpListLength,double* wt);
int getWeight(int **data2,int hapN,int *index,int *snpList,int snpListLength,double *wt, double corrLevel);
int getWeight(int **data2,int hapN,int *index,int *snpList,int snpListLength,double *wt, double corrLevel,int phaseIndicator);

int phase2loci_em(int *vec1,int *vec2,int hapN,double* pAB,double* pA,double* pB);
double logL_em(double **f, int **count);

int readinGenoEigenFormat(FILE *eigenGenoInput,char **genoEigen,int eigenIndN);
int readinSNPEigenFormat(FILE *eigenSNPInput,char **snpNameEigen,int *snpChrEigen,double *snpGPEigen,int *snpPosEigen,char **snpAlleleEigen);
int readinIndEigenFormat(FILE *eigenIndInput,char **indIDEigen,char *indSexEigen,char **indPopEigen);
int readinGenoXpclrFormat(FILE *xpclrGenoInput,char **genoXpclr, int xpclrIndN);
int readinSNPXpclrFormat(FILE *xpclrSNPInput, char *snpNameXpclr,int *snpChrXpclr,double *snpGPXpclr,int *snpPPXpclr,char**snpAlleleXpclr);
double logPSelection4ss_power(int **data1,int **data2,int hapN1,int hapN2,double w,double *pos,int *n1, int mutpos,double *maxS,double *sValue,int sLength,int *snpList,int SnpListLength,double *wt,int *snpv1_1,int *snpv2_1,int *snpv1_2,int *snpv2_2,int *indexHC,double *p2,double *p1, int *x1);
double logPNeutrality4ss_power(double w,double *p2, int *x1,int *n1, int *snpList, int SnpListlength, double *wt,int mutPos);
double get_t(double s,double p0,double pt);
double partA_pAB(double t);
double partA_pAb(double t);
double func4ss_pA(double p);
double func4ss_pA_all(double p);
double prob_softSweep_fun(int n1,int x1,double pB0,double pBt,double pAB0,double pAb0,double c,double s,double w);
void get_haplotypeFreq(double *pB2,double *pAB2,double *pAb2,double *pB1,double *pAB1,double *pAb1,int *snpv1_1,int *snpv2_1,int hapN1,int *snpv1_2,int *snpv2_2,int hapN2,int *x1_n);

#endif //__p2sweeplib__
