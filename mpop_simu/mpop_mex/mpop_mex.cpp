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
gsl_rng *rx;



void mexFunction( int nlhs, mxArray *plhs[], 
		          int nrhs, const mxArray*prhs[] )
     
{ 
    mwSize m1,n1,m2,n2;
    mxClassID  category1,category2,category3;
        
    /* Check for proper number of arguments */
    if (nrhs != 3) { 
	mexErrMsgTxt("Three input arguments required."); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments."); 
    }
	
	/*
	if (! mxIsStruct (prhs[1])){
		mexErrMsgTxt ("expects struct");
	}
	*/

    /* Check  type*/ 
	category1 = mxGetClassID(prhs[0]);
	category2 = mxGetClassID(prhs[1]);
	category3 = mxGetClassID(prhs[2]);
    if (category1!=mxLOGICAL_CLASS){
	mexErrMsgTxt("First input logical is required."); 
    } else if (category2!=mxDOUBLE_CLASS) {
	mexErrMsgTxt("Second input double is required."); 
    }else if (category3!=mxSTRUCT_CLASS) {
	mexErrMsgTxt("Third input structrue is required."); 
    }


    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    m1 = mxGetM(prhs[0]); 
    n1 = mxGetN(prhs[0]);
    m2 = mxGetM(prhs[1]);
    n2 = mxGetN(prhs[1]);
    mexPrintf("%d %d\n",n1,n2);
	if (n1!=n2){
		mexErrMsgTxt("Equal lengths required.");
	}

  T = gsl_rng_ranlxd2;
  rx = gsl_rng_alloc (T);

  int fidx=mxGetFieldNumber(prhs[2],"g");
  if (fidx<0){ ngen=0;}else{
	  ngen = (int)*mxGetPr(mxGetFieldByNumber(prhs[2], 0, fidx));
  }
  mexPrintf("g=%d\n",ngen);

  fidx=mxGetFieldNumber(prhs[2],"N");
  if (fidx<0){ N=0;}else{
	N = (int)*mxGetPr(mxGetFieldByNumber(prhs[2], 0, fidx));
  }
  mexPrintf("N=%d\n",N);

  fidx=mxGetFieldNumber(prhs[2],"e");
  if (fidx<0){ }else{
	seed = (int)*mxGetPr(mxGetFieldByNumber(prhs[2], 0, fidx));
  }
  mexPrintf("e=%d\n",seed);


  if (seed > 0){
    gsl_rng_set(rx,(int)seed);
  } else{
    seed = (int)time(0);
    gsl_rng_set(rx, seed);
  }
  mexPrintf("seed=%d\n",seed);


  Population pop(prhs[0],prhs[1],prhs[2]);
  mexPrintf("N=%d pop.size=%d\n",N,pop.size());

  if( N != pop.size()){
     Population tmppop(pop);
     bool still_there = false;
     while(still_there == false && pop.selected_startfreq>0){
        tmppop = pop;
        still_there = tmppop.resize(rx, N);
     }
     if (!still_there){
      tmppop.resize(rx, N);
     }
     pop = tmppop;
  }
  mexPrintf("N=%d pop.size=%d\n",N,pop.size());
  select_flag = true;
  pop.add_selected_site(rx);
  bool fixed = false;
  Population evolvepop(N);
  evolvepop = pop;
  while(fixed == false && ngen>0){
   evolvepop = pop;
   fixed = evolvepop.evolve(rx, ngen);
  }
  evolvepop.printms(outfile, seed);


   set<double> sites;
   sites = evolvepop.get_sites();
   int count = sites.size();
    double *yp; 
    plhs[1] = mxCreateDoubleMatrix(1, count, mxREAL); 
    yp = mxGetPr(plhs[1]);
    
	bool *yp2;
 	plhs[0] = mxCreateLogicalMatrix(pop.size(), count);
	//yp2 = (bool *)mxGetPr(plhs[0]);
	yp2=(bool *)mxGetLogicals(plhs[0]);
    evolvepop.printms(yp2);

    
	set<double>::iterator it;
    count=0;
	for (it = sites.begin(); it != sites.end(); it++){
		double forout = *it;
		yp[count]=forout;
		count++;
    }

  return;
}

/*
int mainx(int argc, char *argv[]){
  //    get arguments  [msout,N,m,r,s,h,g]

  int i;
  for (i=1; i<argc; i++){
    if(argv[i][0] == '-'){
      switch(argv[i][1]){
      case 'S':
	select_flag = true;
	break;
      case 'O':
	s_p = atof(argv[i]+3);
	break;
      case 's':
	s = atof(argv[i]+3);
	change_s = true;
	break;
      case 'h':
	h = atof(argv[i]+3);
	change_h = true;
	break;
      case 'N':
	N = atoi(argv[i]+3);
	change_N = true;
	break;
      case 'g':
	ngen = atoi(argv[i]+3);
	break;
      case 'r':
	rrate = atof(argv[i]+3);
	change_r = true;
	break;
      case 'm':
	mu =  atof(argv[i]+3);
	change_mu = true;
	break;
      case 'i':
	strcpy(infile, argv[i]+3);
	break;
      case 'o':
	strcpy(outfile, argv[i]+3);
	break;
      case 'e':
       seed = atoi(argv[i]+3);
	break;
      case 'v':
	printf("MPopSim v0.5\n");
	exit(0);
      default:
	printf("Unknown option %s\n", argv[i]);
      }
    }
  }
  T = gsl_rng_ranlxd2;
  rx = gsl_rng_alloc (T);

  if (seed > 0){
    gsl_rng_set(rx,(int)seed);
  }
  else{
    seed = (int)time(0);
    gsl_rng_set(rx, seed);
  }
  Population pop(infile);
  if (change_mu){
    pop.mu = mu;
  }
  if (change_r){
    pop.r = rrate;
  }
  if (change_s){
    pop.set_s(s);
  }
  if (change_h){
    pop.set_h(h);
  }
  if(change_N && N != pop.size()){
    Population tmppop(pop);
    bool still_there = false;
    while(still_there == false && pop.selected_startfreq>0){
      tmppop = pop;
      still_there = tmppop.resize(rx, N);
    }
    if (!still_there){
      tmppop.resize(rx, N);
    }
    pop = tmppop;
  }
  if (select_flag == true){
    pop.add_selected_site(rx);
  }
  if (s_p > 0){
    pop.add_selected_site(s_p);
  }
  bool fixed = false;
  Population evolvepop(N);
  evolvepop = pop;
  while(fixed == false && ngen>0){
   evolvepop = pop;
   fixed = evolvepop.evolve(rx, ngen);
  }
  evolvepop.printms(outfile, seed);
}
  
*/

  
