#ifndef _MONSTER_null_H_
#define _MONSTER_null_H_

#include <math.h>
#include <Eigen>
#include <stdlib.h>
#include <stdio.h>
#include "MatVec.h"
#include "optimize.h"

using namespace Eigen;

typedef Map<Matrix<double, Dynamic, Dynamic, RowMajor> > MapRowXd;
typedef Map<Matrix<double, Dynamic, Dynamic> > MapColXd;
typedef Map<const Matrix<double, Dynamic, Dynamic, RowMajor> > MapRowXdConst;


struct NullResult {
	int n_indiv;
	double *resid; //resid = solve(Sigma_a) * (y-y_hat) / sigma_e
	double **Sigma_inv_half; //Sigma_inv_half = t(Sigma_a)^(-1/2) - (Sigma_a)^(-1) * X * [t(X)*Sigma_a^(-1)*X]^(-1) * t(X) * t(Sigma_a)^(-1/2)
	double **cov;
	double error[3][2];
	double *resid0; //untransformed residual vector
	double *resid_std; //transformed residual vector
};



/////////////////////////////////////////////////////////////////////
/*
	Parameters not to be optimized over in the likelihood function
	U is eigenvectors of phi. phi = U D Ut
*/
struct INFO {
	int n_indiv;
	int n_cov;
	const double *lambda; // eigenvalues of lambda
	const double * U; //eigenvectors of phi, phi = U D Ut. aligned by column of the n_indiv * n_indiv matrix
	double * ytU;
	double *xtU; //aligned by the columns of the n_cov * n_indiv matrix
};
////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////
/*
	Calculations under the null model
*/
///////////////////////////////////////////////////////////////////////////
void SKAT_null (int n_indiv, int n_cov, int res, double *y, double **x, double *lambda, double **U, 
			struct NullResult* nullResult){
			
	double Null_Likelihood (double a, void *info);
	void Null_likelihood_param (double a, void *info, struct NullResult* nullResult, int resid);

	int i;

	/////////////////////////////////////////////////////////
	/// info: extra arguments to feed into Null_likelihood //
    /////////////////////////////////////////////////////////
	struct  INFO *info;	
	info = (struct INFO*) malloc (sizeof(struct INFO));
	info -> n_indiv = n_indiv;
	info -> n_cov = n_cov;	
	info -> lambda = lambda;
	
	// info -> ytU
	MapRowXd U_map (U[0], n_indiv, n_indiv); // Map U to a matrix variables, by rows 
	info -> U = &U_map(0);	// phi = U * diag(lambda) * Ut, U is lined up by the rows
	MapRowXd y_map (y, n_indiv, 1); 
	MatrixXd ytU_map = y_map.transpose() * U_map;
	info -> ytU = &ytU_map(0);

	// info -> xtU, as an 1D array lined up by columns
	MapRowXd x_map(x[0],n_indiv,n_cov); //x_map is a n_indiv * n_cov matrix that shares memory with x	
	MatrixXd xtUmatrix = x_map.transpose() * U_map; //compute xtU
	info -> xtU = &xtUmatrix(0); //xtU is an array that is lined up by the columns of the matrix xtUmatrix

	////////////////////////////////////////////////////////////////
	//// Optimize over exp(-sigma_a2/sigma_e2) within (tol, 1) /////
	////////////////////////////////////////////////////////////////
	double tol =1e-9, optim_a, xmin;	
	optim_a = Brent_fmin(0, 1, Null_Likelihood, info, tol);
	
	////////////////////////////////////////////////////////////////////////////
	// Get quantities under optimal null model. Results stored in nullResult. //
	// See struct NullResult for details                                      //
	////////////////////////////////////////////////////////////////////////////

	Null_likelihood_param (optim_a, info, nullResult, res);

	free(info);
}




////////////////////////////////////////////////////////////////////
/*
	log likelihood function under the null.
	Returns n * log(sigma_e2) + log( det(sigma_a) )
	Input a = exp ( - sigma_a2 / sigma_e2 ) is in (0, 1)
	Requires X to be non-degenerate.
	Phi should be symmetric positive semi-definite?
*/
////////////////////////////////////////////////////////////////////
double Null_Likelihood (double a, void *info){

	a = -log(a); // Now, a = sigma_a2 / sigma_e2
	//read parameters from info
	int n_indiv = ((struct INFO *) info ) -> n_indiv;
	int n_cov = ((struct INFO *) info ) -> n_cov;
	const double *lambda = ((struct INFO *) info ) -> lambda;
	double * ytU = ((struct INFO *) info ) -> ytU;
	double *xtU = ((struct INFO *) info ) -> xtU;	
	MapColXd xtU_map(xtU,n_cov, n_indiv); // xtU_map is n_cov * n_indiv matrix that shares memory with xtU
	
	MapColXd ytU_map(ytU, 1, n_indiv);// map ytU to a matrix variable
	// construct lambda_a_inv
	MatrixXd lambda_a_inv (n_indiv,1); 
	int i;
	for (i=0; i<n_indiv; i++) {lambda_a_inv(i,0) = 1/ (1 + a * lambda[i]);}
	//compute xtSy and sigma_e2
	MatrixXd xtSy = xtU_map * lambda_a_inv.asDiagonal() * ytU_map.transpose();
	MatrixXd xtSx = xtU_map * lambda_a_inv.asDiagonal() *xtU_map.transpose();
	MatrixXd sigma_e2_matrix = ytU_map * lambda_a_inv.asDiagonal() * ytU_map.transpose()-
							xtSy.transpose() * (xtSx).inverse() * xtSy;												
	double sigma_e2 = sigma_e2_matrix(0,0)/n_indiv;	
	// compute likelihood
	for (i=0; i<n_indiv; i++) lambda_a_inv (i, 0) = log(lambda_a_inv(i, 0)); //feed log (lambda_a_inv) into lambda_a_inv
	double likelihood = n_indiv * log(sigma_e2) - lambda_a_inv.sum();

	return (likelihood);
}





/////////////////////////////////////////////////////////////////////////////////
/*
	Compute residuals and Sigma_inv_half under the optimal null model.
	Result is stored in nullResult.
	See definition of struct NullResult for details.
*/
/////////////////////////////////////////////////////////////////////////////////
void Null_likelihood_param (double a, void *info, struct NullResult* nullResult, int res){
	a = -log(a); // Now, a = sigma_a2 / simga_a2
	// Read parameters from info
	int n_indiv = ((struct INFO *) info ) -> n_indiv;
	int n_cov = ((struct INFO *) info ) -> n_cov;
	const double *lambda = ((struct INFO *) info ) -> lambda;
	const double *U = ((struct INFO *) info ) -> U;
	double * ytU = ((struct INFO *) info ) -> ytU;
	double *xtU = ((struct INFO *) info ) -> xtU;

	MapColXd xtU_map(xtU,n_cov, n_indiv); //map xtU to a matrix variable
	MapRowXdConst U_map(U, n_indiv, n_indiv);	// map U to a const matrix varaible, by rows	
	nullResult -> n_indiv = n_indiv; 

	double **Sigma_inv_half, *resid;
	MapColXd ytU_map(ytU, 1, n_indiv); //map ytU to matrix variable
	// construct lambda_a_inv and lambda_a_half
	MatrixXd lambda_a_inv (n_indiv,1);
	MatrixXd lambda_a_half_inv (n_indiv, 1);
	int i, j;
	for (i=0; i<n_indiv; i++) {
		lambda_a_inv(i,0) = 1/ (1 + a * lambda[i]);
		lambda_a_half_inv(i,0) = sqrt (lambda_a_inv(i,0));
	}
	// Compute tmp and tmp2 as matrix variables using Eigen
	// Get Sigma_inv_half_tmp and resid_tmp as matrix variables via Eigen
	// Then compute Sigma_inv_half and resid coarsely using loops
	MatrixXd xtSx = xtU_map * lambda_a_inv.asDiagonal() * xtU_map.transpose();
	MatrixXd tmp = MatrixXd::Identity (n_indiv, n_indiv) - 
				xtU_map.transpose() * ( xtSx ).inverse() * xtU_map * 
				lambda_a_inv.asDiagonal();
	// This is a time consuming step. Takes 24s for n_indiv = 3671
	MatrixXd Sigma_inv_half_tmp = tmp;
	Sigma_inv_half = d2array(n_indiv, n_indiv);
	for (i=0; i<n_indiv; i++){
		for (j=0; j<n_indiv; j++){
			Sigma_inv_half [i][j] = Sigma_inv_half_tmp (i, j) * lambda_a_half_inv (i, 0);
		}
	}
	MatrixXd tmp2= lambda_a_inv.asDiagonal() * tmp * ytU_map.transpose();
	MatrixXd sigma_e2_matrix = ytU_map * tmp2 / n_indiv;
	double sigma_e = sqrt(sigma_e2_matrix(0,0));
	MatrixXd resid_tmp = U_map * tmp2;
	resid = d1array(n_indiv);
	for (i=0; i<n_indiv; i++) resid[i] = resid_tmp(i,0) / sigma_e;
	
	nullResult -> Sigma_inv_half = Sigma_inv_half;
	nullResult -> resid = resid;

	double sigma_a2, sigma_e2, H2, Iaa=0, Iee=0, Iae=0, sd_sigma_a2, sd_sigma_e2, sd_H2;
	sigma_e2 = pow(sigma_e, 2);
	sigma_a2 = sigma_e2 * a;
	H2 = sigma_a2 / (sigma_a2+sigma_e2);
	for (i=0; i<n_indiv; i++){
		Iaa += .5 * pow( lambda[i] / (sigma_a2 * lambda[i] + sigma_e2), 2 );
		Iee += .5 * pow( 1 / (sigma_a2 * lambda[i] + sigma_e2), 2 );
		Iae += .5 * pow( 1 / (sigma_a2 * lambda[i] + sigma_e2), 2 ) * lambda[i];
	}
	double det = Iaa * Iee - pow(Iae, 2);
	double crossterm = -Iae/det;
	sd_sigma_a2 = sqrt(Iee / det);
	sd_sigma_e2 = sqrt(Iae / det);
	(nullResult->error)[1][1] = sd_sigma_a2;//sd_sigma_a2
	(nullResult->error)[2][1] = sd_sigma_e2;//sd_sigma_e2 
	(nullResult->error)[0][1] = sqrt( pow(sigma_e2*sd_sigma_a2, 2) + pow(sigma_a2*sd_sigma_e2, 2) - 2*sigma_e2*sigma_a2*crossterm) / pow(sigma_a2+sigma_e2, 2); //sd_H2
	(nullResult->error)[0][0] = H2;
	(nullResult->error)[1][0] = sigma_a2;
	(nullResult->error)[2][0] = sigma_e2;
	
	nullResult->cov = d2array(n_cov, 2);
	MatrixXd beta = ( xtSx ).inverse() * xtU_map * lambda_a_inv.asDiagonal() * ytU_map.transpose();
	MatrixXd beta_cov = (xtU_map * lambda_a_inv.asDiagonal() * xtU_map.transpose()).inverse();
	for(i=0; i<n_cov; i++){
		(nullResult->cov)[i][0] = beta(i,0);
		(nullResult->cov)[i][1] = sqrt(beta_cov(i,i)) * sigma_e;
	}
	
	if(res==1){
		(nullResult->resid0) = d1array(n_indiv);
		(nullResult->resid_std) = d1array(n_indiv);
		MatrixXd resid_prim(n_indiv, 1);
		MatrixXd resid_std(n_indiv, 1);
		for (i=0; i<n_indiv; i++) {
			resid_prim(i,0) = tmp2(i,0) / lambda_a_inv(i,0);
			resid_std(i,0) = tmp2(i,0) / (lambda_a_half_inv(i,0) * sigma_e);
		}
		resid_prim = U_map * resid_prim;
		resid_std = U_map * resid_std;
		for(i=0; i<n_indiv; i++) {
			(nullResult->resid0)[i] = resid_prim(i, 0);
			(nullResult->resid_std)[i] = resid_std(i,0);
		}
	}	
}


#endif
