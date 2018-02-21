#ifndef _MONSTER_TEST_H_
#define _MONSTER_TEST_H_

#include <math.h>
#include <Eigen>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include "MatVec.h"
#include "MONSTER_null.h"
#include "qfc.h"

using namespace Eigen;

/*
	See the following definitions in "MONSTER_null.h"
		typedef Map<Matrix<double, Dynamic, Dynamic, RowMajor> > MapRowXd;
		typedef Map<Matrix<double, Dynamic, Dynamic> > MapColXd;
		struct NullResult;
		struct INFO;
*/

// Parameters for each Q. used for p values of each Q, p_min and q_min
struct Each_Q_Param {
	double muQ;
	double varQ;
	double df;
};

//
struct INFO_integrate {
	int n_rank; // rank of Z'(I-M)Z (number of nonzeror eigenvalues)
	int n_grid;
	double *rho;
	double *q_min;
	double *lambda;
	double *tau;
	double mu_Q;
	double sigma; // sqrt(sigma_Q^2 - sigma_zeta^2)/sigma_Q
};
	



////////////////////////////////////////////////////////////////////////////////////////////
/*
	Main routine 
	Input: 
	G: double** array, coded by 0/1/2. Each row is an individual, each column is a site.
	W: double* array of leanth n_site.
	nullResult: of class struct NullResult. 
	Return:
	pvalue: the final pvalue as a double.
*/
/////////////////////////////////////////////////////////////////////////////////////////////
double SKAT (int n_indiv, int n_site, int n_grid, double **G, double *W, struct NullResult *nullResult, double *UtGW, double *rho_optim) {
// UtGW is a 1D array storing the n_indiv*n_site matrix "UtGW" by columns
	void Get_Q (int n_grid, int n_indiv, int n_site, double **G, double *W, double *resid, double *rho, double *Q);
	void Each_Q (int n_grid, int n_indiv, int n_site, double *Z, double *rho, double *Q, double *q_min, double *rho_optim);
	void Get_INFO_integrate (int n_indiv, int n_site, int n_grid, double *Z, double *rho, struct INFO_integrate *info_integrate);
	double Get_Pvalue (int n_indiv, int n_site, int n_grid, double *Z, double *rho, double *q_min);

	// Get rho vector
	int i;
	double *rho = d1array(n_grid);
	for (i=0; i<n_grid; i++) {
		rho[i] = (i+1-1.0)/(n_grid-1);
		if (rho[i] > 0.9999) rho[i] = 0.9999;
	}
	
	// Get Z = AGW
	double **Sigma_inv_half = nullResult -> Sigma_inv_half;		
	MapRowXd A_map(Sigma_inv_half[0], n_indiv, n_indiv); // Map Sigma_inv_half to a matrix variable, stored by row.
	MapColXd UtGW_map(UtGW, n_indiv, n_site);
	MatrixXd Z_map = A_map * UtGW_map; // Z_map is a matrix variable, lined up by the columns of Z = AGW.

	// Get the Q's and q_min for each Q	
	double *resid = nullResult -> resid; // Extract resid;
	double *Q = d1array(n_grid);
	Get_Q(n_grid, n_indiv, n_site, G, W, resid, rho, Q);// Compute each Q value
	double *q_min = d1array(n_grid);
	Each_Q (n_grid, n_indiv, n_site, &Z_map(0), rho, Q, q_min, rho_optim); //Get q_min for each Q. &Z_map(0) is a 1D array lined up by the columns of Z
	
	// Get final p value
	double *tmp = &Z_map(0);
	double pvalue = Get_Pvalue (n_indiv, n_site, n_grid, &Z_map(0), rho, q_min);
	
	free(rho);
	free(Q);
	free(q_min);
	return pvalue;

}

///////////////////////////////////////////////////
/*
	Compute each Q value (One for each rho)
*/
///////////////////////////////////////////////////
void Get_Q (int n_grid, int n_indiv, int n_site, double **G, double *W, double *resid, double *rho, double *Q){	
	MapRowXd G_map (G[0], n_indiv, n_site);
	MapRowXd W_map (W, n_site, 1);
	MapRowXd resid_map(resid, n_indiv, 1);
	
	MatrixXd residGW = resid_map.transpose() * G_map * W_map.asDiagonal();
	double tmp1 = pow(residGW.sum(), 2);
	double tmp2 = ( residGW.array() * residGW.array() ).sum();
		
	int i;
	for (i=0; i<n_grid; i++) Q[i] = rho[i] * tmp1 + (1 - rho[i]) * tmp2;	
}

/////////////////////////////////////////////////////////////////////
/*
	Compute the q_min's	(oen for each rho).
	Result stored in double *q_min of length n_grid.
	double *Z is a 1D array lined up by the columns of ""Z = AGW".
*/
/////////////////////////////////////////////////////////////////////
void Each_Q (int n_grid, int n_indiv, int n_site, double *Z, double *rho, double *Q, double *q_min, double *rho_optim){

	void Get_each_Q_Param (int n_grid, int n_indiv, int n_site, double *Z, double *rho, 
			struct Each_Q_Param *each_Q_Param);
	// Obtain for each Q the parameters necessary for comuting the p values and q_min's using Liu's approximation.
	struct Each_Q_Param *each_Q_Param = (struct Each_Q_Param *) malloc (n_grid * sizeof(struct Each_Q_Param)); // To store parameters for each Q
	Get_each_Q_Param (n_grid, n_indiv, n_site, Z, rho, each_Q_Param); 	
	
	// Compute the p values of the Qs using Liu's chi-squared approximation of quadratic forms, matching the 1st four moments.
	int i;
	MatrixXd p_each_Q (n_grid, 1);
	double Q_norm, q_org;
	struct Each_Q_Param *tmp;
	for (i=0; i<n_grid; i++){
		tmp = each_Q_Param + i;
		Q_norm = ( Q[i] - tmp -> muQ ) / sqrt(tmp -> varQ) * sqrt(2 * tmp -> df) + tmp -> df;
		p_each_Q(i,0) = gsl_cdf_chisq_Q (Q_norm, tmp -> df);
	}

	// The smallest p value is p_min.
	MatrixXf::Index mini, minj;
	double p_min = p_each_Q.minCoeff(&mini, &minj);	
	if(rho[mini]>=0.9999){
		(*rho_optim) = 1;
	} else{
		(*rho_optim) = rho[mini];
	}

	// Compute the q_min's, again using Liu's method.
	for (i=0; i<n_grid; i++){
		tmp = each_Q_Param + i;
		q_org = gsl_cdf_chisq_Qinv (p_min, tmp -> df);
		q_min [i] = ( q_org - tmp -> df) / sqrt(2 * tmp -> df) * sqrt(tmp -> varQ) + tmp -> muQ;
	}
	free(each_Q_Param);

}

////////////////////////////////////////////////////////////////////////////////////////////
/*
	Obtain parameters for each Q neccessary for computing the p values and q_min's, using
	Liu's moment matching approximation for quadratic forms.
	Result stored in each_Q_Param (a struct Each_Q_Param * of length n_grid, malloc'ed prior to this routine)/ 
*/
////////////////////////////////////////////////////////////////////////////////////////////
void Get_each_Q_Param (int n_grid, int n_indiv, int n_site, double *Z, double *rho, struct Each_Q_Param *each_Q_Param){
			
	void Get_Liu_Param (double *moments, struct Each_Q_Param *result);
	
	MapColXd Z_map (Z, n_indiv, n_site);
	
	////////////////////////////////////////////////////////////////////////////////////////
	////// P as orthonormal eigenvectors of 1*1t, and  K = t(ZP )* ZP//////
	////////////////////////////////////////////////////////////////////////////////////////	
	
	// Construct P
	MatrixXd P (n_site, n_site);
	int i, j, k;
	for (i=0; i<n_site; i++) P(i,0) = 1 / sqrt(n_site);
	for (i=0; i<n_site; i++){
		for (j=1; j<i+1; j++){
			P(i,j) = -1 / sqrt((n_site - j) * (n_site - j + 1));
		}
		if (i<n_site-1) {
			P(i,i+1) = (n_site - i - 1) / sqrt((n_site - i - 1) * (n_site -i));
		}
		for (j=i+2; j<n_site; j++){
			P(i,j)=0;
		}
	}
	
	//////////////////////////////////////
	//// Obtain parameters for each Q ////
	//////////////////////////////////////
	
	// Compute K = t(ZP )* ZP;
	MatrixXd ZP = Z_map * P;
	MatrixXd K = ZP.transpose() * ZP;

	// get 1st 4 moments
	MatrixXd DKD = K * 1;
	double D1, D2;
	double moments[4];	
	for (i=0; i<n_grid; i++){
		//get DKD matrix
		DKD = K * 1;
		D1 = n_site * rho[i] + 1 - rho[i];
		D2 = 1- rho[i];
		for (j=1; j<n_site; j++){
			DKD(j,0) *= sqrt(D1 * D2);
			DKD(0,j) *= sqrt(D1 * D2);
			for (k=1; k<n_site; k++){
				DKD(j,k) *= D2;
			}
		}
		DKD(0,0) *= D1;
				
		// eigenvalues of DKD
		SelfAdjointEigenSolver<MatrixXd> es(DKD, EigenvaluesOnly);		

		// 1st 4 moments of the distribution of the Q corresponding to rho[i]
		moments[0] = es.eigenvalues().sum();
		moments[1] = ( es.eigenvalues().array().abs() * es.eigenvalues().array().abs()  ).sum();
		moments[2] = ( es.eigenvalues().array().abs()  * es.eigenvalues().array().abs()  * es.eigenvalues().array().abs()  ).sum();
		moments[3] = ( es.eigenvalues().array().abs()  * es.eigenvalues().array().abs()  * es.eigenvalues().array().abs()  * es.eigenvalues().array().abs()  ).sum();		
		// get parameters for the Q corresponding to rho[i]
		// Stored in the "struct Each_Q_Param" pointer "each_Q_param + i"		
		Get_Liu_Param (moments, each_Q_Param + i); 
	}
	
}


/////////////////////////////////////////////////////////////////////////////////////////////
/*
	H.Liu, Y.Tang, and H.H.Zhang. A new chi-square approximation to the distribution of 
	non-negative definite quadratic forms in non-central normal variables.
	From "Function.R" in SKAT package in R.
	Input:
	"moments" is a 1D double array of length 4, containing the 1st four moments of the targeted quadratic form.
	"result" stores the result. It is a struct Each_q_Param* of length 1.
*/
void Get_Liu_Param (double *moments, struct Each_Q_Param *result){
	double muQ = moments[0];
	double sigmaQ = sqrt (2 * moments[1]);
	double s1 = moments[2] / pow(moments[1], 1.5);
	double s2 = moments[3] / pow(moments[1], 2);
	
	double a, d, l;
	if (pow(s1, 2) > s2){
		a = 1 / ( s1 - sqrt( pow(s1, 2) - s2 ) );
		d = s1 * pow(a, 3) - pow(a, 2);
		l = pow(a, 2) - 2 * d;
	} else {
		l = 1 / s2;
	}
	
	result -> muQ = muQ;
	result -> varQ = pow(sigmaQ, 2);
	result -> df = l;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Computes the final p value.
	Input:
	Z is a 1D double array, lined up by the columns of "Z = AGW".
	q_min is 1D double array of length n_grid, containing the p_min quantiles of the null distributions of the Q's
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Get_Pvalue (int n_indiv, int n_site, int n_grid, double *Z, double *rho, double *q_min){

	void Get_INFO_integrate (int n_indiv, int n_site, int n_grid, double *Z, double *rho, struct INFO_integrate *info_integrate);
	double Ff (double x, void *info_integrate);

	// Get parameters (mu_Q, sigma_Q, sigma_zeta, tau_rho ...) necessary as inputs of the Ff function to be integrated
	// See def of struct INFO_integrate for details.
	struct INFO_integrate *info_integrate = (struct INFO_integrate *) malloc (sizeof(struct INFO_integrate));
	Get_INFO_integrate (n_indiv, n_site, n_grid, Z, rho, info_integrate);
	info_integrate -> q_min = q_min;
	
	// Numerical integration of Ff from 0 to +infinity.
	int iter_limit = 10000;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (iter_limit); // Allocate workspace   
    double result, result_new, error, relerr = 1e-3;   
	int status = 1;
    gsl_function F;
    F.function = &Ff;
    F.params = info_integrate;    
	gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
	gsl_integration_qagiu (&F, 0, 0, relerr, iter_limit, w, &result, &error); 
	if (result > 0.9){
		relerr= 1e-4;
		status = gsl_integration_qagiu (&F, 0, 0, relerr, iter_limit, w, &result_new, &error); 							 
		if (status == 0) result = result_new;
	}
	if (result > 0.95){
		status = 1;
		relerr = 1e-6;
		while(status && relerr < 1e-4) {
			status = gsl_integration_qagiu (&F, 0, 0, relerr, iter_limit, w, &result_new, &error); 							 
			relerr *= 5;
        }
		if (status == 0) result = result_new;	
	}
	gsl_set_error_handler(old_handler); 	
	
	double pvalue = 1- result;
	gsl_integration_workspace_free (w);
	
	free(info_integrate -> lambda);
	free(info_integrate -> tau);
	free(info_integrate);	
	return pvalue;

}


/////////////////////////////////////////////////
/*
	Get parameters necessary as inputs of the Ff function to be integrated.
	Result is sotred in struct INFO_integrate *info_integrate (length 1). See details in the def of this structure.
*/
/////////////////////////////////////////////////
void Get_INFO_integrate (int n_indiv, int n_site, int n_grid, double *Z, double *rho, struct INFO_integrate *info_integrate){

	// Compute matrices ZtMZ, ZtI_MZ, ZtZ1, and the quantity ztz.
	int i;
	MapColXd Z_map (Z, n_indiv, n_site);
	MatrixXd ZtZ = Z_map.transpose() * Z_map;
	VectorXd z = Z_map.rowwise().sum();
	MatrixXd ZtZ1 = Z_map.transpose() * z;
	double ztz = (z.array() * z.array()).sum();
	MatrixXd ZtMZ = ZtZ1 * ZtZ1.transpose() / ztz;
	MatrixXd ZtI_MZ = ZtZ - ZtMZ;	

	////////////////////////////////
	//// Get n_rank and lambda /////
	////////////////////////////////
	SelfAdjointEigenSolver<MatrixXd> es(ZtI_MZ, EigenvaluesOnly); // Eigen-decomposition without eigenvectors.
	int n_rank = 0; // The number of positive eigenvalues of ZtI_MZ.
	// Compute n_rank.
	for (i=0; i<n_site; i++){
		if (es.eigenvalues()(i) >= es.eigenvalues().mean() / 100000) {
			n_rank ++;
		} else if (es.eigenvalues()(i) < -pow(10,-5) ) {
//			printf("ERROR: ZtI_MZ has a negative eigenvalue smaller than -1e-5: %lf!\n",es.eigenvalues()(i));
//			exit(1);
		}
	}
	// Copy positive entries of es.eigenvalues() to lambda.
	int flag = 0;
	double *lambda = d1array(n_rank);
	for (i=0; i<n_site; i++){
		if (es.eigenvalues()(i) >= es.eigenvalues().mean() / 100000 ) {
			lambda[flag] = es.eigenvalues()(i);
			flag++;
		}
	}
	
	////////////////////////////
	//// Get mu_Q and sigma ////
	////////////////////////////
	double sigma_zeta2, mu_Q, sigma_Q2, sigma;		
	sigma_zeta2 = 4.0 * (ZtMZ * ZtI_MZ).trace();	
	mu_Q = es.eigenvalues().sum();
	sigma_Q2 = 2.0 * (es.eigenvalues().array() * es.eigenvalues().array()).sum() + sigma_zeta2;
	sigma = sqrt(1.0 - sigma_zeta2 / sigma_Q2); 

	// get tau
	double *tau;
	tau = d1array(n_grid);
	double tmp = (ZtZ1.array() * ZtZ1.array()).sum() / ztz;
	for (i=0; i<n_grid; i++){
		tau[i] = rho[i] * ztz + (1.0 - rho[i]) * tmp;
	}	
	
	info_integrate -> n_rank =n_rank;
	info_integrate -> n_grid = n_grid;
	info_integrate -> rho = rho;
	info_integrate -> lambda = lambda;
	info_integrate -> tau = tau;
	info_integrate -> mu_Q = mu_Q;
	info_integrate -> sigma = sigma;
	
}

//////////////////////////////////////////////////////
/*
	Function to be integrated to get the final pvalue.
	Ff = F( delta(x) | lambda) * f( x | chi-squared).
	x is the argument to integrated over.
	info_integrate contains other parameters.
*/
//////////////////////////////////////////////////////

double Ff (double x, void *info_integrate){

	// Read parameters from info_integrate.
	int n_rank = ((struct INFO_integrate *)info_integrate) -> n_rank;
	int n_grid =  ((struct INFO_integrate *)info_integrate) -> n_grid;
	double *rho = ((struct INFO_integrate *)info_integrate) -> rho;
	double *q_min = ((struct INFO_integrate *)info_integrate) -> q_min;
	double *lambda = ((struct INFO_integrate *)info_integrate) -> lambda;
	double *tau = ((struct INFO_integrate *)info_integrate) -> tau;
	double mu_Q = ((struct INFO_integrate *)info_integrate) -> mu_Q;
	double sigma = ((struct INFO_integrate *)info_integrate) -> sigma;	
	MapRowXd rho_map (rho, n_grid, 1);
	MapRowXd q_min_map (q_min, n_grid, 1);
	MapRowXd tau_map (tau, n_grid, 1);
	
	// f
	double f = gsl_ran_chisq_pdf (x, 1.0);
	
	// F
	double delta = ( ((q_min_map.array() - tau_map.array() * x) / (1-rho_map.array())).minCoeff() - mu_Q ) * sigma + mu_Q; // delta(x)
	// F ( delta(x) | lambda) using qfc
	double *noncentral = d1array(n_rank);
	int *df = i1array(n_rank);
	int i;
	for (i=0; i<n_rank; i++) df[i] = 1;
	double s =0; // sigma in qfc
	int lim = 100000;
	double acc = 0.000001;
	double trace[7];
	int ifault;
	double F = -1;
	qfc (lambda, noncentral, df, &n_rank, &s, &delta, &lim, &acc, trace, &ifault, &F);
	// Force F to be 0 or 1 for small or large values. Otherwise, qfc will return 0.5 for +/-infinity or extreme positive values.
	// Also force F to be 0 if "q_min_map(n_grid-1) - tau_map(n_grid-1) * x < 0", because in that case delta is -infinity
	if(delta<0) {
		F = 0;
	} else if(delta > pow(10,80)){
		F = 1;
	}
	if (F<0) F=0;
	if(F>1) F=1;
	
	free(noncentral);
	free(df);
	return (F*f);
}



double burdenTest(int n_indiv, int n_site, double **G, double *W, struct NullResult *nullResult, double *UtGW){
	double **Sigma_inv_half = nullResult -> Sigma_inv_half;		
	double *resid = nullResult -> resid;
	MapRowXd G_map (G[0], n_indiv, n_site);
	MapRowXd W_map (W, n_site, 1);
	MapRowXd A_map (Sigma_inv_half[0], n_indiv, n_indiv);
	MapRowXd resid_map (resid, n_indiv, 1);
	MapColXd UtGW_map(UtGW, n_indiv, n_site);
	
	MatrixXd residGW = resid_map.transpose() * G_map * W_map.asDiagonal();
	double Q = pow(residGW.sum(), 2);

	MatrixXd Z_map = A_map * UtGW_map; // Z_map is a matrix variable, lined up by the columns of Z = AGW.
	double coef = (Z_map.rowwise().sum().array() * Z_map.rowwise().sum().array()).sum();
	double pvalue = gsl_cdf_chisq_Q (Q/coef, 1);
	return pvalue;
}


double SecondOrderTest(int n_indiv, int n_site, double **G, double *W, struct NullResult *nullResult, double *UtGW){
	double **Sigma_inv_half = nullResult -> Sigma_inv_half;		
	double *resid = nullResult -> resid;
	MapRowXd G_map (G[0], n_indiv, n_site);
	MapRowXd W_map (W, n_site, 1);
	MapRowXd A_map (Sigma_inv_half[0], n_indiv, n_indiv);
	MapRowXd resid_map (resid, n_indiv, 1);
	MapColXd UtGW_map(UtGW, n_indiv, n_site);

	MatrixXd residGW = resid_map.transpose() * G_map * W_map.asDiagonal();
	double Q = ( residGW.array() * residGW.array() ).sum();

	MatrixXd Z_map = A_map * UtGW_map; // Z_map is a matrix variable, lined up by the columns of Z = AGW.
	MatrixXd ZtZ = Z_map.transpose() * Z_map;
	SelfAdjointEigenSolver<MatrixXd> es(ZtZ, EigenvaluesOnly);	
	
	int i;
	double *lambda = d1array(n_site);
	double *noncentral = d1array(n_site);
	int *df = i1array(n_site);	
	for (i=0; i<n_site; i++){
		lambda[i] = es.eigenvalues()(i);
		df[i] = 1;
	}	
	double s =0; // sigma in qfc
	int lim = 1000000;
	double acc = 0.00000001;
	double trace[7];
	int ifault;
	double pvalue = -1;
	qfc (lambda, noncentral, df, &n_site, &s, &Q, &lim, &acc, trace, &ifault, &pvalue);
	pvalue = 1- pvalue;
	
	free(lambda);
	return pvalue;
}


#endif