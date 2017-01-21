#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include "population.h"
#include <list>
#include <math.h>
using namespace std;
#include "mex.h"

/* Input Arguments */

#define	MS_IN	prhs[0]
#define	PR_IN	prhs[1]

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
    double *yp; 
    double *t,*y; 
    mwSize m,n;
    mxClassID  category1,category2;
    


    
    /* Check for proper number of arguments */
    
    if (nrhs != 2) { 
	mexErrMsgTxt("Two input arguments required."); 
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 

    /* Check  type*/ 
	category1 = mxGetClassID(prhs[0]);
	category2 = mxGetClassID(prhs[1]);
    if (category1~=mxLOGICAL_CLASS){
	mexErrMsgTxt("First input logical is required."); 
    } else if (category2 ~=mxSTRUCT_CLASS) {
	mexErrMsgTxt("Second input structrue is required."); 
    } 


    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) || 
	(MAX(m,n) != 4) || (MIN(m,n) != 1)) { 
	mexErrMsgTxt("YPRIME requires that Y be a 4 x 1 vector."); 
    } 

    /* set parameters */ 	
	// analyze_structure(prhs[1]);


  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc (T);

  if (seed > 0){
    gsl_rng_set(r,(int)seed);
  }
  else{
    seed = (int)time(0);
    gsl_rng_set(r, seed);
  }
  Population pop(prhs[0],prhs[1]);

  bool fixed = false;
  Population evolvepop(N);
  evolvepop = pop;
  while(fixed == false && ngen>0){
   evolvepop = pop;
   fixed = evolvepop.evolve(r, ngen);
  }
  evolvepop.printms(outfile, seed);




    /* Create a matrix for the return argument */ 
    YP_OUT = mxCreateDoubleMatrix(m, n, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    yp = mxGetPr(YP_OUT);
    
    t = mxGetPr(T_IN); 
    y = mxGetPr(Y_IN);
        
    /* Do the actual computations in a subroutine */
    yprime(yp,t,y);
    return;    
}







static void
analyze_structure(const mxArray *structure_array_ptr)
{
  mwSize total_num_of_elements;
  mwIndex index;
  int number_of_fields, field_index;
  const char  *field_name;
  const mxArray *field_array_ptr;
  

  mexPrintf("\n");
  total_num_of_elements = mxGetNumberOfElements(structure_array_ptr); 
  number_of_fields = mxGetNumberOfFields(structure_array_ptr);
  
  /* Walk through each structure element. */
  for (index=0; index<total_num_of_elements; index++)  {
    
    /* For the given index, walk through each field. */ 
    for (field_index=0; field_index<number_of_fields; field_index++)  {
      mexPrintf("\n\t\t");
      display_subscript(structure_array_ptr, index);
         field_name = mxGetFieldNameByNumber(structure_array_ptr, 
                                             field_index);
      mexPrintf(".%s\n", field_name);
      field_array_ptr = mxGetFieldByNumber(structure_array_ptr, 
					   index, 
					   field_index);
      if (field_array_ptr == NULL) {
	mexPrintf("\tEmpty Field\n");
      } else {
         /* Display a top banner. */
         mexPrintf("------------------------------------------------\n");
         get_characteristics(field_array_ptr);
         analyze_class(field_array_ptr);
         mexPrintf("\n");
      }
    }
      mexPrintf("\n\n");
  }
  
  
}




/* -------------------------------- */

int main(int argc, char *argv[]){
  /* 
     get arguments
  */

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
  r = gsl_rng_alloc (T);

  if (seed > 0){
    gsl_rng_set(r,(int)seed);
  }
  else{
    seed = (int)time(0);
    gsl_rng_set(r, seed);
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
      still_there = tmppop.resize(r, N);
    }
    if (!still_there){
      tmppop.resize(r, N);
    }
    pop = tmppop;
  }
  if (select_flag == true){
    pop.add_selected_site(r);
  }
  if (s_p > 0){
    pop.add_selected_site(s_p);
  }
  bool fixed = false;
  Population evolvepop(N);
  evolvepop = pop;
  while(fixed == false && ngen>0){
   evolvepop = pop;
   fixed = evolvepop.evolve(r, ngen);
  }
  evolvepop.printms(outfile, seed);
}
  
  
  
