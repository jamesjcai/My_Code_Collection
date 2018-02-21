

###############################
# Beta function

Generate_Covariates<-function(n){
	
	X1<-rnorm(n)
	X2<-rbinom(n,1,0.5)
	
	return(cbind(X1,X2))

}

Get_Alpha0_C<-function(X){
	
	n<-dim(X)[1]
	re<- X[,1] * 0.5 / sqrt(2) + Beta0 /2
	return(re)
}	

Get_Alpha0_Q<-function(X){

	re<- X[,1] * 0.5 + X[,2] * 0.5 
	return(re)
}	



Get_BETA_C<-function( IDX, MAF_Rare_c, MAF_Common_c){


	n.MAF_Rare_c <-length(MAF_Rare_c)
	n.MAF_Common_c<-length(MAF_Common_c)

	beta_Rare<-rep(0,n.MAF_Rare_c)
	beta_Common<-rep(0,n.MAF_Common_c)
		
	if(IDX == 1 || IDX==9){
		beta_Rare[1:n.MAF_Rare_c]<-log(2)
		beta_Common[1:n.MAF_Common_c]<-log(1.1)	
	} else if(IDX == 2){
		beta_Rare[1:n.MAF_Rare_c]<-log(2)
		beta_Common[1:n.MAF_Common_c]<-log(1.2)	
	} else if (IDX==3){

		beta_Rare[1:n.MAF_Rare_c]<- - 0.2 * log10(MAF_Rare_c)
		beta_Common[1:n.MAF_Common_c]<-- 0.2 * log10(MAF_Common_c)	
		
	} else if (IDX==4){

		beta_Rare[1:n.MAF_Rare_c]<- - 0.4 * log10(MAF_Rare_c)
		beta_Common[1:n.MAF_Common_c]<-- 0.4 * log10(MAF_Common_c)	
		
	} else if (IDX == 5){
		beta_Rare[1:n.MAF_Rare_c]<-log(2)
		
	} else if (IDX == 6){

		beta_Common[1:n.MAF_Common_c]<-log(1.2)	
	} else if (IDX == 7){
	
		beta_Rare[1:n.MAF_Rare_c]<- - 0.4 * log10(MAF_Rare_c)
		beta_Common[1:n.MAF_Common_c]<-- 0.4 * log10(MAF_Common_c)	
		
		p.neg<-round((n.MAF_Rare_c + n.MAF_Common_c) * 0.3)
		idx<-sample(1:(n.MAF_Rare_c + n.MAF_Common_c), p.neg)
		
		id1<-which(idx <= n.MAF_Rare_c)
		id2<-which(idx > n.MAF_Rare_c)
		
		if(length(id1) > 0){
		
			id1<-idx[id1]
			beta_Rare[id1] = beta_Rare[id1] * (-1)
		}
		
		if(length(id2) > 0){
			id2<-idx[id2] -n.MAF_Rare_c 
			beta_Common[id2] = beta_Common[id2] * (-1)
		}
	
	} else if (IDX == 8){
	
		beta_Rare[1:n.MAF_Rare_c]<-log(2)
		beta_Common[1:n.MAF_Common_c]<-log(1.2)	
		
		p.neg<-round((n.MAF_Rare_c + n.MAF_Common_c) * 0.3)
		idx<-sample(1:(n.MAF_Rare_c + n.MAF_Common_c), p.neg)
		
		id1<-which(idx <= n.MAF_Rare_c)
		id2<-which(idx > n.MAF_Rare_c)
		
		if(length(id1) > 0){
		
			id1<-idx[id1]
			beta_Rare[id1] = beta_Rare[id1] * (-1)
		}
		
		if(length(id2) > 0){
			id2<-idx[id2] -n.MAF_Rare_c 
			beta_Common[id2] = beta_Common[id2] * (-1)
		}
	
	}
	
	
	return(list(beta_Rare=beta_Rare, beta_Common=beta_Common))

}

Get_Causal<-function( causal_p, MAF_Rare, MAF_Common, idx1){

	#causal_p<-Causal_percent
	MAF_A<-c(MAF_Rare, MAF_Common)
	n.MAF <-length(MAF_A)
	n.Rare<-length(MAF_Rare)
	n.Common<-length(MAF_Common)
		
	n1<-round(n.MAF * causal_p)
	
	IDX_C<-sort(sample(1:n.MAF, n1))
	IDX_Rare_C1<-which(IDX_C <= n.Rare)
	IDX_Common_C1<-which(IDX_C > n.Rare) 
	
	if(length(IDX_Common_C1) == 0){
		
		IDX_C<-c(IDX_C, sample((n.Rare+1):n.MAF, 1) )
		IDX_Common_C1<-which(IDX_C > n.Rare) 
		
	}

	
	IDX_Rare_C<-IDX_C[IDX_Rare_C1]
	IDX_Common_C<-IDX_C[IDX_Common_C1] - n.Rare
	
	if(idx1== 9){
		causal_p.c<-0.5
		n1<-round(n.Rare * causal_p)
		n2<-round(n.Common * causal_p.c)
		
		if(n1==0){
			n1=1
		}
		if(n2==0){
			n2=1
		}
		
		IDX_Rare_C<-sort(sample(1:n.Rare, n1))
		IDX_Common_C<-sort(sample(1:n.Common, n2))
	}


	MAF_Rare_C<-MAF_Rare[IDX_Rare_C]
	MAF_Common_C<-MAF_Common[IDX_Common_C]
	
	return(list(IDX_Rare_C = IDX_Rare_C, IDX_Common_C = IDX_Common_C, MAF_Rare_C = MAF_Rare_C, MAF_Common_C = MAF_Common_C))
	
}

