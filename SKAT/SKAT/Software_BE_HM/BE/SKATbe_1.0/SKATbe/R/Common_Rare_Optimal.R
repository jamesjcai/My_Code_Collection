
#
#	Function get parameters of optimal test
#
SKAT_2Kernel_Ortho_Optimal_Param<-function(Z1, Z2, r.all){


	n<-dim(Z1)[1]
	p.m<-dim(Z1)[2]
	r.n<-length(r.all)
	
	lambda<-list()
	
	lambda[[1]]<-SKAT:::Get_Lambda(t(Z1) %*% Z1)
	lambda[[2]]<-SKAT:::Get_Lambda(t(Z2) %*% Z2)
	
	par.moments<-list()
	c1<-rep(0,4)
	for(i in 1:2){
		lambda.temp<-lambda[[i]] 
		c1[1]<-sum(lambda.temp)
		c1[2]<-sum(lambda.temp^2)
		c1[3]<-sum(lambda.temp^3)
		c1[4]<-sum(lambda.temp^4)
		param.temp<-SKAT:::Get_Liu_Params_Mod(c1)

		muQ<-param.temp$muQ
		varQ<-param.temp$sigmaQ^2
		df<-param.temp$l
		
		par.moments[[i]]<-c(muQ, varQ, df)

	}

	out<-list(lambda = lambda, par.moments=par.moments)
	
	return(out)
}



#
#	Function get SKAT statistics with given rho
#		Q.all is a matrix with n.q x n.r

SKAT_Optiaml_Each_Q<-function(param.m, Q.all, r.all, c1.all){

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]

	pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
	pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
	param.mat<-NULL

	for(i in 1:n.r){
		Q<-Q.all[,i]
		r.corr<-r.all[i]
		
		c1<-c1.all[,i]
		param.temp<-SKAT:::Get_Liu_Params_Mod(c1)

		muQ<-param.temp$muQ
		varQ<-param.temp$sigmaQ^2
		df<-param.temp$l

		# get pvalue
		Q.Norm<-(Q - muQ)/sqrt(varQ) * sqrt(2*df) + df
		pval[,i]<- 1-pchisq(Q.Norm,  df = df)

		param.mat<-rbind(param.mat,c(muQ,varQ,df))
		
		#cat(c(i, Q, muQ, varQ, df))
		#cat("\n")
		#cat(lambda.temp)
		#cat("\n")
	}
	
	pmin<-apply(pval,1,min)
	for(i in 1:n.r){
	
		muQ<-param.mat[i,1]
		varQ<-param.mat[i,2]
		df<-param.mat[i,3]

		q.org<-qchisq(1-pmin,df=df)
		q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
		pmin.q[,i]<-q.q

	}
	
	# pmin : min p-value
	# pmin.q : q-values of min p.
	out<-list(pmin=pmin,pval=pval,pmin.q=pmin.q)
	return(out)

}


#
#	x is the common varinat term
# x is q2
SKAT_2Kernel_Ortho_Optimal_Integrate_Func<-function(x,pmin.q,param.m,r.all){

	par1<-param.m$par.moments[[1]]
	par2<-param.m$par.moments[[2]]
	
	x1<-(x - par2[3])/ sqrt(2* par2[3]) * sqrt(par2[2]) + par2[1] 
	n.r<-length(r.all)
	n.x<-length(x)
	
	temp1<-r.all %x% t(x1)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)
	
	

		

	temp.q1<-(temp.min - par1[1])/sqrt(par1[2])*sqrt(2*par1[3]) + par1[3]
	#temp.q2<-(x - par2[1])/sqrt(par2[2]) * sqrt(2* par2[3]) + par2[3]
	temp.q2<-x
	re<-pchisq(temp.q1 ,df=par1[3]) * dchisq(temp.q2,df=par2[3])
	#re<-dchisq(temp.q2, df=par2[3])
	
	x<<-x
	r.all<<-r.all
	temp<<-temp
	temp1<<-temp1
	temp.min<<-temp.min
	pmin.q<<-pmin.q
	
	temp.q1<<-temp.q1
	temp.q2<<-temp.q2
	re<<-re
	
	#cat("ZZZ", x ,"\n")
	#cat("AAA", temp.min ,"\n")
	#cat("BBB", c(temp.q1), "\n")
	#cat("CCC", c(temp.q2), "\n")
	#cat("DDD", max(re), min(x), max(x) , "\n")
	return(re)
	

}



SKAT_2Kernel_Ortho_Optimal_PValue_Liu<-function(pmin.q, param.m, r.all){

	#print(pmin.q)
	#print(param.m)
	#cat("Param end \n")
	#cat("Upper=", upper, " Lower=", lower,"\n")
	
	par2<-param.m$par.moments[[2]]
	df = par2[3]
	upper<<- 30 * df
	 re<-integrate(SKAT_2Kernel_Ortho_Optimal_Integrate_Func, lower=0, upper=upper, subdivisions=500
	,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15)
	
	pvalue<-1-re[[1]]
	return(pvalue)

}




SKAT_2Kernel_Optimal_Get_Q<-function(Z1, Z2, res, r.all, n.Resampling, res.out, res.moments=NULL){

	n.r<-length(r.all)
	p.m<-dim(Z1)[2]

	Q.r<-rep(0,n.r)
	Q.r.res<-NULL
	Q.sim<-NULL	
	
	temp1<-t(res) %*% Z1
	temp2<-t(res) %*% Z2
	for(i in 1:n.r){
		r.corr<-r.all[i]
		Q1<-(1-r.corr) * rowSums(temp1^2)
		Q2<-r.corr * rowSums(temp2^2)
		Q.r[i]<-Q1 + Q2
	}
	Q.r = Q.r /2
  	if(n.Resampling > 0){
	
		temp1<-t(res.out) %*% Z1
		temp2<-t(res.out) %*% Z2
		Q.r.res<-matrix(rep(0,n.Resampling *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-r.all[i]
			Q1<-(1-r.corr) * rowSums(temp1^2)
			Q2<-r.corr * rowSums(temp2^2)
			Q.r.res[,i]<-Q1 + Q2
		}
		Q.r.res = Q.r.res/2
  	}

	if(!is.null(res.moments)){

		temp1<-t(res.moments) %*% Z1
		temp2<-t(res.moments) %*% Z2
		n.moments<-dim(res.moments)[2]
		Q.sim<-matrix(rep(0,n.moments *n.r),ncol=n.r)
		for(i in 1:n.r){
			r.corr<-r.all[i]
			Q1<-(1-r.corr) * rowSums(temp1^2)
			Q2<-r.corr * rowSums(temp2^2)
			Q.sim[,i]<-Q1 + Q2
		}
		Q.sim = Q.sim/2

	}

	re<-list(Q.r=Q.r, Q.r.res=Q.r.res , Q.sim=Q.sim)
	return(re)

}

######################################
# Get Params with a given r
#
#	(A+B)^2 = AA + AB + BA + BB 
#	(A+B)^3 = AAA + ABA + BAA + BBA + AAB + ABB + BAB + BBB
#	(A+B)^4 = AAAA + ABAA + BAAA + BBAA + AABA + ABBA + BABA + BBBA 
#			+ AAAB + ABAB + BAAB + BBAB + AABB + ABBB + BABB + BBBB
#
#	tr[(A+B)^2] = AA + 2AB + BB
#	tr[(A+B)^3] = AAA + 3 AAB + 3 BBA + BBB
#	tr[(A+B)^4] = AAAA + BBBB + 4 AAAB + 4 BBBA + 4 AABB + 2 ABAB 

SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r<-function(Z1.1, Z2.1, r.all){

	c1<-matrix(rep(0,4* length(r.all)), ncol=length(r.all))
	
	A1<-t(Z1.1) %*% Z1.1
	B1<-t(Z2.1) %*% Z2.1
	
	A2<-A1 %*% A1
	B2<-B1 %*% B1
	
	A11<-t(Z1.1) %*% Z2.1
	A22<-A11 %*% t(A11)
	B22<-t(A11) %*% A11
	B333<-t(A11) %*%  A1 %*% A11
	
	#####################################
	#
	
	c1[1,]<-sum(Z1.1^2) * (1-r.all) + sum(Z2.1^2) * r.all
	c1[2,]<-sum(A1^2) * (1-r.all)^2 + sum(B1^2) * (r.all)^2 + sum(A11^2) * 2 * (1-r.all) * r.all
	c1[3,]<-sum(A2 * A1) * (1-r.all)^3 + sum(B2 * B1) * (r.all)^3 + sum(A22 * A1) * 3 * (1-r.all)^2 * r.all + sum(B1 * B22) * 3 * (1-r.all) * r.all^2
	c1[4,]<-sum(A2 * A2) * (1-r.all)^4 + sum(B2 * B2) * (r.all)^4 + sum(A22 * A2) * 4 * (1-r.all)^3 * r.all + sum(B2 * B22) * 4 * (1-r.all) * r.all^3 + sum(B1 * B333) * 4 * (1-r.all)^2 * r.all^2 + sum(B22 * B22) * 2 * (1-r.all)^2 * r.all^2
				
	
	#c2<-matrix(rep(0,4* length(r.all)), ncol=length(r.all))
	
	#for(i in 1:4){
	#	for(j in 1:length(r.all)){
	#		c2[i,j]<-sum(lambda.all[[j]]^i)
	#	}
	#
	#}
	
	#cat(c1)
	#cat(c2)
	#cat(sum((c1-c2)^2))
	
	#c1<<-c1
	#c2<<-c2
	return(c1)

}


SKAT_2Kernel_Ortho_Optimal_Get_Pvalue<-function(Q.all, Z1.1, Z2.1, r.all, method){

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]
	p.m<-dim(Z1)[2]

	
	#lambda.all<-list()
	#
	# Not computationally effecient. 
	# Should be changed!!
	
	#K1<-Z1.1 %*% t(Z1.1)
	#K2<-Z2.1 %*% t(Z2.1)
	#for(i in 1:n.r){
	#	r.corr<-r.all[i]
	#	K<- (1-r.corr) * K1 + r.corr * K2
	#	lambda.all[[i]]<-SKAT:::Get_Lambda(K)
	#}
	#c1.all<-matrix(rep(0,4* length(r.all)), ncol=length(r.all))
	#for(i in 1:4){
	#	for(j in 1:length(r.all)){
	#		c1.all[i,j]<-sum(lambda.all[[j]]^i)
	#	}
	#
	#}

	c1.all<-SKAT_2Kernel_Ortho_Optimal_Get_Params_each_r(Z1.1, Z2.1, r.all)


	# Get Mixture param 
	param.m<-SKAT_2Kernel_Ortho_Optimal_Param(Z1.1, Z2.1,r.all)
	Each_Info<-SKAT_Optiaml_Each_Q(param.m, Q.all, r.all, c1.all)
	pmin.q<-Each_Info$pmin.q
	pval<-rep(0,n.q)


	for(i in 1:n.q){
		pval[i]<-SKAT_2Kernel_Ortho_Optimal_PValue_Liu(pmin.q[i,], param.m, r.all)
	}
		
	return(list(p.value=pval,p.val.each=Each_Info$pval))

}

#######################################################33
#	Linear

# Z1: Genotypes of rare variants
# Z2: Genotypes of common variants

SKAT_CR_Optimal = function(res,Z1, Z2, X1, kernel=NULL, weights = NULL, s2 = 0, pi_1=NULL ,method=NULL
, res.out=NULL, n.Resampling =0, r.all, out_type="C" , burden1 = FALSE, burden2 = FALSE){

	# if r.all >=0.999 ,then r.all = 0.999. It is just for computation.
	IDX<-which(r.all >= 0.999)
	if(length(IDX) > 0){
		r.all[IDX]<-0.999
	}

	n<-dim(Z1)[1]
	p.m1<-dim(Z1)[2]
	p.m2<-dim(Z2)[2]	
	n.r<-length(r.all)

	# Weights for rare variants, but no weights for common variants	
	if (kernel == "linear.weighted") {
		Z1 = t(t(Z1) * (weights))
	}
	
 	if(out_type =="C"){
	
		Z1.1<-Z1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z1)
		Z2.1<-Z2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z2)
	
	} else {
	
		Z1.1<- (Z1 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z1 * pi_1))
		Z2.1<- (Z2 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z2 * pi_1))
	}


	temp1<-t(Z1.1) %*% Z1.1
	temp2<-t(Z2.1) %*% Z2.1
	z1.w<-sum(diag(temp1))
	z2.w<-sum(diag(temp2))
	
	
	#Z1<-Z1/sqrt(z1.w) * sqrt(z1.w + z2.w)
	#Z2<-Z2/sqrt(z2.w) * sqrt(z1.w + z2.w)

	z1.var<-sum(temp1 * temp1)
	z2.var<-sum(temp2 * temp2)
	
	
	Z1<-Z1/sqrt(z1.var)
	Z2<-Z2/sqrt(z2.var)

	# Make Orthogonal
	out.QR<-qr(Z1)
	Z1.Q<-try(as.matrix(qr.Q(out.QR)), silent = TRUE)
	Z1.R<-try(as.matrix(qr.R(out.QR)), silent = TRUE)
	
	if(class(Z1.Q) == "try-error" || class(Z1.R) == "try_error"){

		cat("A")
		out.QR<-qr(Z1, LAPACK = TRUE)	
		
		Z1.Q<-try(as.matrix(qr.Q(out.QR)), silent = TRUE)
		Z1.R<-try(as.matrix(qr.R(out.QR)), silent = TRUE)

		if(class(Z1.Q) == "try-error" || class(Z1.R) == "try_error"){

			stop("QR decomposition error!")	
		
		}
	}

	
	Z2.item<-t(Z1.Q) %*% Z2
	Z2.item1<- Z1.Q %*% Z2.item
	Z2.Ortho<-Z2 - Z2.item1
	
	Z2<-Z2.Ortho

	if(burden1){
		Z1<-as.matrix(rowSums(Z1))
	}
	if(burden2){
		Z2<-as.matrix(rowSums(Z2))
	}

  	#Get P0Z, wher P0 = diag(n) - X1%*%solve( t(X1)%*%X1)%*%t(X1)
  	
  	if(out_type =="C"){
	
		Z1.1<-Z1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z1)
		Z2.1<-Z2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z2)
	
	} else {
	
		Z1.1<- (Z1 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z1 * pi_1))
		Z2.1<- (Z2 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z2 * pi_1))
		s2<-1
	}
	
	###########################################
	# Compute Q.r and Q.r.res
	##########################################
	out.Q<-SKAT_2Kernel_Optimal_Get_Q(Z1, Z2, res, r.all, n.Resampling, res.out)
	Q.all<-rbind(out.Q$Q.r, out.Q$Q.r.res) / s2

	##################################################
	# Compute P-values 
	#################################################

	out<-SKAT_2Kernel_Ortho_Optimal_Get_Pvalue(Q.all, Z1.1 / sqrt(2), Z2.1 / sqrt(2), r.all, method)
	
	param<-list(p.val.each=NULL,q.val.each=NULL)
	param$p.val.each<-out$p.val.each[1,]
	param$q.val.each<-Q.all[1,]
	param$rho<-r.all
	param$minp<-min(param$p.val.each)

	id_temp<-which(param$p.val.each == min(param$p.val.each))
	id_temp1<-which(param$rho >= 0.999) # treat rho > 0.999 as 1
	if(length(id_temp1) > 0){
		param$rho[id_temp1] = 1
	}
	param$rho_est<-param$rho[id_temp]


	p.value<-out$p.value[1]

	p.value.resampling<-NULL
	if(n.Resampling > 0){
		p.value.resampling<-out$p.value[-1]
		param$pval.each.resample<-out$p.val.each[-1,]
	}

 	re<-list(p.value = p.value, p.value.resampling = p.value.resampling
	, Test.Type = method, Q = NA, param=param )  
  	
	return(re)	

}



SKAT_CR_Linear = function(obj.res,Z1, Z2, kernel="linear.weighted", weights.beta=c(1,25), weights.beta2=c(1,1),
weights = NULL, r.all=c(0,0.25,0.5,0.75,1), burden1 = FALSE, burden2 = FALSE){

	if(is.null(weights)){
		MAF<-colMeans(Z1)/2
		weights<-SKAT:::Beta.Weights(MAF,weights.beta)
	}
	MAF2<-colMeans(Z2)/2
	weights2<-SKAT:::Beta.Weights(MAF2,weights.beta2)
	Z2 = t(t(Z2) * (weights2))
	
	if(obj.res$out_type == "C"){
		s2 = obj.res$s2
		pi_1 = NULL
		out_type="C" 
	} else {
	
		s2 = 1
		pi_1 = obj.res$pi_1	
		out_type="D"
	}
	
	
	re<-SKAT_CR_Optimal(obj.res$res,Z1, Z2, obj.res$X1, kernel, weights, s2=s2, pi_1=pi_1, method=NULL,
	res.out= obj.res$res.out, n.Resampling = obj.res$n.Resampling, r.all = r.all, out_type=out_type, 
	burden1=burden1, burden2=burden2)

	return(re)
}


SKAT_CR_Linear_Reverse = function(obj.res,Z1, Z2, kernel="linear.weighted", weights.beta=c(1,25), weights.beta2=c(1,1), 
weights = NULL, r.all=c(0,0.25,0.5,0.75,1), burden1 = FALSE, burden2 = FALSE){

	if(is.null(weights)){
		MAF<-colMeans(Z1)/2
		weights<-SKAT:::Beta.Weights(MAF,weights.beta)
		if (kernel == "linear.weighted") {
			Z1 = t(t(Z1) * (weights))
		}
	
	}
	MAF2<-colMeans(Z2)/2
	weights2<-SKAT:::Beta.Weights(MAF2,weights.beta2)
	Z2 = t(t(Z2) * (weights2))
	
	if(obj.res$out_type == "C"){
		s2 = obj.res$s2
		pi_1 = NULL
		out_type="C" 
	} else {
	
		s2 = 1
		pi_1 = obj.res$pi_1	
		out_type="D"
	}
	
	
	re<-SKAT_CR_Optimal(obj.res$res,Z2, Z1, obj.res$X1, kernel="linear", weights, s2=s2, pi_1=pi_1, method=NULL,
	res.out= obj.res$res.out, n.Resampling = obj.res$n.Resampling, r.all = r.all, out_type=out_type, 
	burden1=burden1, burden2=burden2)

	return(re)
}



SKAT_Scale_Genotypes= function(obj.res, Z1, Z2, weights.beta=c(1,25), weights.beta2=c(1,1), weights=NULL, r.corr1, r.corr2){

	if(obj.res$out_type == "C"){
		s2 = obj.res$s2
		pi_1 = NULL
		out_type="C" 
	} else {
	
		s2 = 1
		pi_1 = obj.res$pi_1	
		out_type="D"
	}
	
	X1=obj.res$X1

	n<-dim(Z1)[1]
	p.m1<-dim(Z1)[2]
	p.m2<-dim(Z2)[2]	

	if(is.null(weights)){
		MAF<-colMeans(Z1)/2
		weights<-SKAT:::Beta.Weights(MAF,weights.beta)
	}
	
	# Weights for rare variants, but no weights for common variants	
	Z1 = t(t(Z1) * (weights))
	
	MAF2<-colMeans(Z2)/2
	weights2<-SKAT:::Beta.Weights(MAF2,weights.beta2)
	Z2 = t(t(Z2) * (weights2))
	
	
   # r.corr
   if(r.corr1 == 1){
  		Z1<-cbind(rowSums(Z1))
   } else if(r.corr1 > 0){

   		p.m<-dim(Z1)[2]	
		R.M<-diag(rep(1-r.corr1,p.m)) + matrix(rep(r.corr1,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z1<- Z1 %*% t(L) 
   }

   if(r.corr2 == 1){
  		Z2<-cbind(rowSums(Z2))
   } else if(r.corr2 > 0){

   		p.m<-dim(Z2)[2]	
		R.M<-diag(rep(1-r.corr2,p.m)) + matrix(rep(r.corr2,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z2<- Z2 %*% t(L) 
   }
  
  	if(out_type =="C"){
	
		Z1.1<-Z1 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z1)
		Z2.1<-Z2 - X1%*%solve( t(X1)%*%X1)%*%(t(X1) %*% Z2)
	
	} else {
	
		Z1.1<- (Z1 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z1 * pi_1))
		Z2.1<- (Z2 * sqrt(pi_1)) - (X1 * sqrt(pi_1))%*%solve(t(X1)%*%(X1 * pi_1))%*% (t(X1) %*% (Z2 * pi_1))
		s2<-1
	}
  
	temp1<-t(Z1.1) %*% Z1.1
	temp2<-t(Z2.1) %*% Z2.1
	z1.mean<-sum(diag(temp1))
	z2.mean<-sum(diag(temp2))
	
	z1.var<-sum(temp1 * temp1)
	z2.var<-sum(temp2 * temp2)
	
	
	
	Z1<-Z1/sqrt(sqrt(z1.var))
	Z2<-Z2/sqrt(sqrt(z2.var))
	
	return(list(Z1=Z1, Z2=Z2, z1.mean=z1.mean, z2.mean=z2.mean, z1.var=z1.var, z2.var=z2.var))
}



SKAT_Joint_CommonRare<-function(Z, obj, weights.beta.rare=c(1,25), weights.beta.common=c(1,1), 
r.corr.rare=0, r.corr.common=0, Cutoff = NULL, is_dosage=FALSE, SetID=NULL){


	#
	if(is.matrix(Z) != TRUE){
		stop("Z should be a matrix")
	}

	# Compute common and rare
	n <- dim(Z)[1]
    m <- dim(Z)[2]
    
    if(is.null(Cutoff)){
    	Cutoff<-1/sqrt(N.Sample * 2)
    }
    # Check Cutoff
    if(Cutoff < 0 && Cutoff > 0.5){
    	stop("Error in Cutoff! It should be NULL or a numeric value between 0 to 0.5")
    }
    
	
	MAF<-SKAT:::Get_MAF(Z)
	id.rare<-intersect(which(MAF < Cutoff), which(MAF > 0))
	id.common<-intersect(which(MAF >= Cutoff), which(MAF > 0))
	
	n.rare= length(id.rare)
	n.common= length(id.common)
			
	if(n.rare == 0 && n.common > 0){
		
		# Run SKAT with common variants
		re<-SKAT_With_NullModel(Z, obj, weights.beta=weights.beta.common, weights = NULL, r.corr=r.corr.common
		, is_check_genotype=TRUE, is_dosage = is_dosage, missing_cutoff=1, SetID = SetID)

	} else if (n.common == 0 && n.rare > 0){
		
		# Run SKAT with rare variants
		re<-SKAT_With_NullModel(Z, obj, weights.beta=weights.beta.rare, weights = NULL, r.corr=r.corr.rare
		, is_check_genotype=TRUE, is_dosage = is_dosage, missing_cutoff=1, SetID = SetID)
		
	} else {
		
		Z.rare<-cbind(Z[,id.rare])
		Z.common<-cbind(Z[,id.common])
		
		obj1<-SKAT_Scale_Genotypes(obj, Z.rare, Z.common, r.corr1=r.corr.rare, r.corr2=r.corr.common, 
		weights.beta=weights.beta.rare, weights.beta2=weights.beta.common )
		
		re<-SKAT(cbind(obj1$Z1, obj2$Z2), kernel="linear", obj.res)
		
	}
	
	re$n.rare = n.rare
	re$n.common = n.common
	re$method= "Joint.CommonRare"
	
	return(re)

}
