/*=================================================================
 *
 * YPRIME.C	Sample .MEX file corresponding to YPRIME.M
 *	        Solves simple 3 body orbit problem 
 *
 * The calling syntax is:
 *
 *		[yp] = yprime(t, y)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2006 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.10.6.4 $ */
#include <math.h>
#include "mex.h"

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

static	double	mu = 1/82.45;
static	double	mus = 1 - 1/82.45;


static void yprime(
		   double	yp[],
		   double	*t,
 		   double	y[]
		   )
{
    double	r1,r2;
    
    (void) t;     /* unused parameter */

    r1 = sqrt((y[0]+mu)*(y[0]+mu) + y[2]*y[2]); 
    r2 = sqrt((y[0]-mus)*(y[0]-mus) + y[2]*y[2]);

    /* Print warning if dividing by zero. */    
    if (r1 == 0.0 || r2 == 0.0 ){
	mexWarnMsgTxt("Division by zero!\n");
    }
    
    yp[0] = y[1];
    yp[1] = 2*y[3]+y[0]-mus*(y[0]+mu)/(r1*r1*r1)-mu*(y[0]-mus)/(r2*r2*r2);
    yp[2] = y[3];
    yp[3] = -2*y[1] + y[2] - mus*y[2]/(r1*r1*r1) - mu*y[2]/(r2*r2*r2);
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		          int nrhs, const mxArray*prhs[] )
     
{ 
    double *yp; 
    double *t,*y; 
    mwSize m,n; 
mwSize i;
    
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    m = mxGetM(prhs[0]); 
    n = mxGetN(prhs[0]);
    
    //times=new double[m];

    /* Create a matrix for the return argument */ 
    YP_OUT = mxCreateDoubleMatrix(m, n, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    yp = mxGetPr(YP_OUT);
    
    t = mxGetPr(prhs[0]); 
    y = mxGetPr(Y_IN);

	//memcpy(t,times,m*sizeof(double));
    //memcpy(mxGetPr(prhs[0]), y_pr, m * sizeof(double));
	for (i=0; i<n; i++){
		mexPrintf("%f\n", t[i]);
	}
    
    /* Do the actual computations in a subroutine */
    // yprime(yp,t,y); 
    return;
    
}


