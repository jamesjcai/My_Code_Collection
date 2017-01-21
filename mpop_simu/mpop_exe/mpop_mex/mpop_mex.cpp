#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include "population.h"
#include <list>
#include <math.h>
#include "mex.h"
using namespace std;

/* Input Arguments */

#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]


/* Output Arguments */

#define	YP_OUT	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/*
  mpop v0.5: simulates a constant size population under a Wright-Fisher 
  model

  Options:
    -g number of generations
    -N number of chromosomes
    -r recombination rate from one end of the locus to the other
    -m mutation rate locus-wide
    -i input file (either a single ms population or an output files from mpop)
    -S flag whether or not a new selected mutation is added (always added at pos 0.5)
    -O adds an old mutation to the population at a given frequency
    -s
    -h s and h are the selection coefficients on the selected locus  w11= 1, w12 = 1+hs, w22 = 1+s
    -e random number seed

    TO DO:
    -allow for multiple selected sites?
    -model migration?
*/


int N=0; 
bool change_N = false;
int ngen=0;
double rrate;
bool change_r = false; 
double mu; 
bool change_mu = false;
double s;
bool change_s = false; 
double h;
bool change_h = false;
bool select_flag = false;
char infile[50]; 
char outfile[50] = "mpopout";
int seed = 0;
double s_p = 0;
const gsl_rng_type *T;
gsl_rng *r;



void mexFunction( int nlhs, mxArray *plhs[], 
    	     	  int nrhs, const mxArray*prhs[] )
     
{ 

/*void mexFunction(int nlhs,mxArray *[],int nrhs,const mxArray *prhs[])*/

    /* Check for proper number of input and output arguments */    
    if (nrhs != 1) {
	mexErrMsgTxt("One input argument required.");
    } 
    if (nlhs > 1){
	mexErrMsgTxt("Too many output arguments.");
    }

}
  
  
  
