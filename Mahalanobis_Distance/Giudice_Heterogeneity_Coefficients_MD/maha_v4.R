

# Two simple R functions that compute Mahalanobis' D, confidence intervals, overlap coefficents, and heterogeneity coefficients, either from raw data with maha() or from Cohen's d values and a common correlation matrix with maha.summary(). For general information see Del Giudice (2009, 2013).

# By Marco Del Giudice (2016). Version 4.

# Note on confidence intervals: exact confidence intervals are computed with Reiser's (2001) method. Especially for small D values, the equations may not be solvable; in those cases, one or both CI bounds are set to NA. For more information see Reiser (2001). Bootstrapped CIs are bias-corrected and accelerated; for details see Kelley (2005).

# Note on heterogeneity coefficients: the heterogeneity coefficient H quantifies heterogeneity in the variables' contribution to the multivariate effect size (0 = max homogeneity; 1 = max heterogeneity). The EPV coefficient (equivalent proportion of variables) is the proportion of contributing variables that would result in the same heterogeneity, in a hypothetical scenario where a certain proportion of variables contribute equally to D while the remaining ones make no contribution. See Del Giudice (in press - 2016).

# References:

# Del Giudice, M. (2009). On the real magnitude of psychological sex differences. Evolutionary Psychology, 7, 264-279. doi:10.1177/147470490900700209
# Del Giudice, M. (2013). Multivariate misgivings: Is D a valid measure of group and sex differences? Evolutionary Psychology, 11, 1067-1076. doi:10.1177/147470491301100511
# Del Giudice, M. (in press - 2016). Heterogeneity coefficients for Mahalanobis' D as a multivariate effect size. Multivariate Behavioral Research.
# Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean difference: Bootstrap and parametric confidence intervals. Educational and Psychological Measurement, 65, 51-69. doi:10.1177/0013164404264850
# Reiser, B. (2001). Confidence intervals for the Mahalanobis distance. Communications in Statistics: Simulation and Computation, 30, 37–45. doi:10.1081/SAC-100001856




# function maha(data_A, data_B, alpha=NULL, conf.level=.95, boot.n=10000)
#
# Returns the Mahalanobis distance D, confidence intervals, heterogeneity coefficients, and coefficients of overlap, computed from raw data. Can compute disattenuated estimates if desired. Note: The correlation matrices of the two groups are pooled by taking weighted averages before computing D. Tucker’s Congruence Coefficient provides an index of similarity between the two correlation matrices (0 = min similarity; 1 = max similarity). 
#
# Arguments
# data_A, data_B: raw matrices/data frames for the two groups. Must contain the same variables.
# alpha: vector of reliability coefficients (optional: only required for disattenuation)
# conf.level: CI width (optional)
# boot.n: number of bootstrap samples (optional; deafult is 10000)
#
# Value
# returns a list object containing some or all of the following:
# D: Mahalanobis' D
# CI_lower_exact: exact CI lower bound (NA if not solvable)
# CI_upper_exact: exact CI upper bound (NA if not solvable)
# CI_lower_boot: bootstrapped CI lower bound 
# CI_upper_boot: bootstrapped CI upper bound 
# OVL: coefficient of overlap (single distribution)
# OVL2: Cohen's coefficient of overlap 1-U1 (joint distribution)
# CC_cor: Tucker's Congruence Coefficient (similarity between correlation matrices, 0-1)
# H: heterogeneity coefficient H (0 = max homogeneity; 1 = max heterogeneity)
# EPV: EPV coefficient (equivalent proportion of variables, 0-1)
# Dc: disattenuated D
# OVLc: disattenuated coefficient of overlap (single distribution)
# OVL2c: disattenuated Cohen's coefficient of overlap (joint distribution)
# Hc: heterogeneity coefficient for disattenuated D 
# EPVc: EPV coefficient for disattenuated D
# d_values: vector of Cohen's d values
# dc_values: vector of disattenuated Cohen's d values




# function maha.summary(d_values, cor_matrix, alpha=NULL, nA=NULL, nB=NULL, conf.level=.95)
#
# Returns the Mahalanobis distance D, confidence intervals, heterogeneity coefficients, and coefficients of overlap, computed from summary statistics. Can compute disattenuated estimates if desired.
#
# Arguments
# d_values: row vector of standardized differences (Cohen's d)
# cor_matrix: correlation matrix
# alpha: vector of reliability coefficients (optional: only required for disattenuation)
# nA, nB: sample size of the two groups (optional: only required for exact CI)
# conf.level: CI width (optional)
#
# Value
# returns a list object containing some or all of the following:
# D: Mahalanobis' D
# CI_lower_exact: exact CI lower bound (NA if not solvable)
# CI_upper_exact: exact CI upper bound (NA if not solvable)
# OVL: coefficient of overlap (single distribution)
# OVL2: Cohen's coefficient of overlap 1-U1 (joint distribution)
# H: heterogeneity coefficient H (0 = max homogeneity; 1 = max heterogeneity)
# EPV: EPV coefficient (equivalent proportion of variables, 0-1)
# Dc: disattenuated D
# OVLc: disattenuated coefficient of overlap (single distribution)
# OVL2c: disattenuated Cohen's coefficient of overlap (joint distribution)
# Hc: heterogeneity coefficient for disattenuated D 
# EPVc: EPV coefficient for disattenuated D 
# d_values: vector of Cohen's d values
# dc_values: vector of disattenuated Cohen's d values




maha <- function(data_A, data_B, alpha=NULL, conf.level=.95, boot.n=10000) {

# preliminary computations: Ns, pooled variances, pooled correlation matrix

nA = as.numeric(length(complete.cases(data_A)))
nB = as.numeric(length(complete.cases(data_B)))

cor_A = cor(data_A, use="pairwise.complete.obs")
cor_B = cor(data_B, use="pairwise.complete.obs")
pooled_cor = (cor_A*nA+cor_B*nB)/(nA+nB)

pooled_variances = ((nA-1)*diag(var(data_A, na.rm=TRUE))+(nB-1)*diag(var(data_B, na.rm=TRUE)))/(nA+nB)

# compute Cohen's d values and Mahalanobis D

d_values = (colMeans(data_A, na.rm=TRUE)-colMeans(data_B, na.rm=TRUE))/sqrt(pooled_variances)

D2 = mahalanobis(d_values, array(0, dim=c(length(d_values))), pooled_cor)
output = list()
output$D = sqrt(D2)
		
# compute exact CI

p = length(d_values)	

Fcal = (D2)*(nA*nB*(nA+nB-p-1))/(p*(nA+nB)*(nA+nB-2))

lower_prob = conf.level+(1-conf.level)/2
upper_prob = (1-conf.level)/2
critical_F_upper = qf(upper_prob, p, nA+nB-p-1)
critical_F_lower = qf(lower_prob, p, nA+nB-p-1)

# lower
if (D2>critical_F_lower) { ncp_est_max=10000/(1/nA+1/nB)
ncp_est_min = 0
ncp_est = (ncp_est_max-ncp_est_min)/2

est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)

while(abs(est_p-lower_prob) > .00001) {
	if ((est_p-lower_prob) < 0) {ncp_est_max=ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	else {ncp_est_min = ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
	}
output$CI_lower_exact = sqrt(ncp_est*(1/nA+1/nB))
	}
else output$CI_lower_exact=NA

# upper
if (D2>critical_F_upper) { ncp_est_max=10000/(1/nA+1/nB)
ncp_est_min = 0
ncp_est = (ncp_est_max-ncp_est_min)/2

est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)

while(abs(est_p-upper_prob) > .00001) {
	if ((est_p-upper_prob) < 0) {ncp_est_max = ncp_est
		ncp_est=ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	else {ncp_est_min = ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
	}
output$CI_upper_exact = sqrt(ncp_est*(1/nA+1/nB))
	}
else output$CI_upper_exact = NA

# compute bootstrapped CI

boot_D = numeric(boot.n)
alpha_level = 1-conf.level

for (sample in 1:boot.n) {
	data_A.sampled = data_A[sample(seq(1:nrow(data_A)), replace = TRUE),]
	data_B.sampled = data_B[sample(seq(1:nrow(data_B)), replace = TRUE),]
	
	nA.sampled = length(complete.cases(data_A.sampled))
	nB.sampled = length(complete.cases(data_B.sampled))
	cor_A.sampled = cor(data_A.sampled, use="pairwise.complete.obs")
	cor_B.sampled = cor(data_B.sampled, use="pairwise.complete.obs")
	pooled_cor.sampled = (cor_A.sampled*nA.sampled+cor_B.sampled*nB.sampled)/(nA.sampled+nB.sampled)
	pooled_variances.sampled = ((nA.sampled-1)*diag(var(data_A.sampled, na.rm=TRUE))+(nB.sampled-1)*diag(var(data_B.sampled, na.rm=TRUE)))/(nA.sampled+nB.sampled)
	d_values.sampled = (colMeans(data_A.sampled, na.rm=TRUE)-colMeans(data_B.sampled, na.rm=TRUE))/sqrt(pooled_variances.sampled)
	boot_D[sample] = sqrt(mahalanobis(d_values.sampled, array(0, dim=c(length(d_values.sampled))), pooled_cor.sampled))
}

jackknife_results = numeric(nrow(data_A)+nrow(data_B))
n.1 = nrow(data_A)
n.2 = nrow(data_B)

Marker.1 = seq(1, n.1, 1)

for (sample in 1:n.1) {
	data_A.jack = data_A[Marker.1[-sample],]
	data_B.jack = data_B
	
	nA.jack = length(complete.cases(data_A.jack))
	nB.jack = length(complete.cases(data_B.jack))
	cor_A.jack = cor(data_A.jack, use="pairwise.complete.obs")
	cor_B.jack = cor(data_B.jack, use="pairwise.complete.obs")
	pooled_cor.jack = (cor_A.jack*nA.jack+cor_B.jack*nB.jack)/(nA.jack+nB.jack)
	pooled_variances.jack = ((nA.jack-1)*diag(var(data_A.jack, na.rm=TRUE))+(nB.jack-1)*diag(var(data_B.jack, na.rm=TRUE)))/(nA.jack+nB.jack)
	d_values.jack = (colMeans(data_A.jack, na.rm=TRUE)-colMeans(data_B.jack, na.rm=TRUE))/sqrt(pooled_variances.jack)
	jackknife_results[sample] = sqrt(mahalanobis(d_values.jack, array(0, dim=c(length(d_values.jack))), pooled_cor.jack))
}

Marker.2 = seq(1, n.2, 1)

for (sample in 1:n.2) {
	data_A.jack = data_A
	data_B.jack = data_B[Marker.2[-sample],]
	
	nA.jack = length(complete.cases(data_A.jack))
	nB.jack = length(complete.cases(data_B.jack))
	cor_A.jack = cor(data_A.jack, use="pairwise.complete.obs")
	cor_B.jack = cor(data_B.jack, use="pairwise.complete.obs")
	pooled_cor.jack = (cor_A.jack*nA.jack+cor_B.jack*nB.jack)/(nA.jack+nB.jack)
	pooled_variances.jack = ((nA.jack-1)*diag(var(data_A.jack, na.rm=TRUE))+(nB.jack-1)*diag(var(data_B.jack, na.rm=TRUE)))/(nA.jack+nB.jack)
	d_values.jack = (colMeans(data_A.jack, na.rm=T)-colMeans(data_B.jack, na.rm=TRUE))/sqrt(pooled_variances.jack)
	jackknife_results[n.1+sample] = sqrt(mahalanobis(d_values.jack, array(0, dim=c(length(d_values.jack))), pooled_cor.jack))
}

Mean.Jackknife = mean(jackknife_results)

a = (sum((Mean.Jackknife-jackknife_results)^3))/(6*sum((Mean.Jackknife-jackknife_results)^2)^(3/2))

z0 = qnorm(sum(boot_D < output$D)/boot.n)

CI.Low.BCa = pnorm(z0 + (z0+qnorm(alpha_level/2))/(1-a*(z0+qnorm(alpha_level/2))))
CI.Up.BCa = pnorm(z0 + (z0+qnorm(1-alpha_level/2))/(1-a*(z0+qnorm(1-alpha_level/2))))

output$CI_lower_boot = quantile(boot_D, CI.Low.BCa, names=FALSE)
output$CI_upper_boot = quantile(boot_D, CI.Up.BCa, names=FALSE)

# compute overlap coefficients

output$OVL = 2*pnorm(-output$D/2)
output$OVL2 = output$OVL/(2-output$OVL)

# compute Tucker's CC

cor_ab = cor_A*cor_B
cor_a2 = cor_A^2
cor_b2 = cor_B^2

output$CC_cor = sum(cor_ab[lower.tri(cor_ab)])/sqrt(sum(cor_a2[lower.tri(cor_a2)])*sum(cor_b2[lower.tri(cor_b2)]))

# compute heterogeneity coefficients

C_values = sort((t(d_values)%*%solve(pooled_cor))*d_values)
C_values[C_values<0] = 0
N = length(C_values)
sum_Ci = sum(C_values)
sum_iCi = sum(C_values*seq(1:N))

output$H = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
output$EPV = 1-(output$H*(N-1)/N)

output$d_values=d_values

# disattenuation

if (is.null(alpha)==FALSE) {
	d_values.d = d_values/sqrt(alpha)
	alpha_matrix = sqrt(crossprod(t(alpha), t(alpha)))
	diag(alpha_matrix) = 1
	pooled_cor.d = pooled_cor/alpha_matrix
	D2.d = mahalanobis(d_values.d, array(0, dim=c(length(d_values.d))), pooled_cor.d)
	output$Dc = sqrt(D2.d)
	
	output$OVLc = 2*pnorm(-output$Dc/2)
	output$OVL2c = output$OVLc/(2-output$OVLc)
	
	C_values = sort((t(d_values.d)%*%solve(pooled_cor.d))*d_values.d)
	C_values[C_values<0] = 0
	N = length(C_values)
	sum_Ci = sum(C_values)
	sum_iCi = sum(C_values*seq(1:N))

	output$Hc = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
	output$EPVc = 1-(output$Hc*(N-1)/N)
	}
	
if (is.null(alpha)==FALSE) {output$dc_values = d_values.d
	}

# output

return(output)

}




maha.summary <- function(d_values, cor_matrix, alpha=NULL, nA=NULL, nB=NULL, conf.level=.95) {

# compute Mahalanobis D

D2 = mahalanobis(d_values, array(0, dim=c(length(d_values))), cor_matrix)
output = list()
output$D = sqrt(D2)		

# compute exact CI

if (is.null(nA)==FALSE)	{p = length(d_values)
	Fcal = (D2)*(nA*nB*(nA+nB-p-1))/(p*(nA+nB)*(nA+nB-2))
	lower_prob = conf.level+(1-conf.level)/2
	upper_prob = (1-conf.level)/2
	critical_F_upper = qf(upper_prob, p, nA+nB-p-1)
	critical_F_lower = qf(lower_prob, p, nA+nB-p-1)

# lower
if (D2>critical_F_lower) { ncp_est_max = 10000/(1/nA+1/nB)
ncp_est_min = 0
ncp_est = (ncp_est_max-ncp_est_min)/2

est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)

while(abs(est_p-lower_prob) > .00001) {
	if ((est_p-lower_prob) < 0) {ncp_est_max = ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	else {ncp_est_min = ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
	}
output$CI_lower_exact = sqrt(ncp_est*(1/nA+1/nB))
	}
else output$CI_lower_exact = NA

# upper
if (D2>critical_F_upper) { ncp_est_max = 10000/(1/nA+1/nB)
ncp_est_min = 0
ncp_est = (ncp_est_max-ncp_est_min)/2

est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est)

while(abs(est_p-upper_prob) > .00001 ) {
	if ((est_p-upper_prob) < 0) {ncp_est_max = ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	else {ncp_est_min = ncp_est
		ncp_est = ncp_est_min+(ncp_est_max-ncp_est_min)/2
		}
	est_p = pf(Fcal, p, (nA+nB-p-1), lower.tail=TRUE, ncp=ncp_est);
	}
output$CI_upper_exact = sqrt(ncp_est*(1/nA+1/nB))
	}
else output$CI_upper_exact = NA
	}

# compute overlap coefficents

output$OVL = 2*pnorm(-output$D/2)
output$OVL2 = output$OVL/(2-output$OVL)

# compute heterogeneity coefficients

C_values = sort((t(d_values)%*%solve(cor_matrix))*d_values)
C_values[C_values<0] = 0
N = length(C_values)
sum_Ci = sum(C_values)
sum_iCi = sum(C_values*seq(1:N))

output$H = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
output$EPV = 1-(output$H*(N-1)/N)

output$d_values = d_values

# disattenuation

if (is.null(alpha)==FALSE) {
	d_values.d = d_values/sqrt(alpha)
	alpha_matrix = sqrt(crossprod(t(alpha), t(alpha)))
	diag(alpha_matrix) = 1
	cor_matrix.d = cor_matrix/alpha_matrix
	D2.d = mahalanobis(d_values.d, array(0, dim=c(length(d_values.d))), cor_matrix.d)
	output$Dc = sqrt(D2.d)

	output$OVLc = 2*pnorm(-output$Dc/2)
	output$OVL2c = output$OVLc/(2-output$OVLc)
	
	C_values = sort((t(d_values.d)%*%solve(cor_matrix.d))*d_values.d)
		C_values[C_values<0] = 0
	N = length(C_values)
	sum_Ci = sum(C_values)
	sum_iCi = sum(C_values*seq(1:N))

	output$Hc = ( (2/N)*sum_iCi -((N+1)/N)*sum_Ci) / ((N-1)*mean(C_values))
	output$EPVc = 1-(output$Hc*(N-1)/N)
	}
	
if (is.null(alpha)==FALSE) {output$dc_values = d_values.d
	}
	
# output

return(output)

}
